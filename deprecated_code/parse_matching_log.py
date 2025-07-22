import sys, re, csv
from pathlib import Path

# ── 1.  CLI & file check ─────────────────────────────────────────────
if len(sys.argv) != 3:
    sys.exit("usage: python parse_matching_log.py <logfile> <output.csv>")

log_path = Path(sys.argv[1])
if not log_path.is_file():
    sys.exit(f"file not found: {log_path}")

out = Path(sys.argv[2])


# ── 2.  Regexes (number first for fallback / no‑match) ───────────────
re_target   = re.compile(r"Target:\s+(\S+)")
re_matched  = re.compile(r"Matched\s+(\d+)/\d+\s+fragments")
re_fallback = re.compile(r"(\d+)\s+matched only with fallback", re.I)
re_nomatch  = re.compile(r"(\d+)\s+had no population fragment", re.I)

summary, current = {}, None
stats = {"matched": 0, "fallback": 0, "no_match": 0}

with log_path.open() as fh:
    for line in fh:
        if m := re_target.search(line):
            if current is not None:
                summary[current] = stats
                stats = {"matched": 0, "fallback": 0, "no_match": 0}
            current = m.group(1)

        elif m := re_matched.search(line):
            stats["matched"] = int(m.group(1))

        elif m := re_fallback.search(line):
            stats["fallback"] = int(m.group(1))

        elif m := re_nomatch.search(line):
            stats["no_match"] = int(m.group(1))

if current is not None:                 # flush last target
    summary[current] = stats

# ── 3.  Write CSV ────────────────────────────────────────────────────
# out = Path("fragment_summary.csv")
with out.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["target", "matched", "fallback", "no_match"])
    for tgt, s in summary.items():
        w.writerow([tgt, s["matched"], s["fallback"], s["no_match"]])

# ── 4.  Console preview ──────────────────────────────────────────────
print(f"📄  Results for {log_path.name}")
print("target          matched  fallback  no_match")
print("──────────────  ───────  ────────  ────────")
for tgt, s in summary.items():
    print(f"{tgt:14}  {s['matched']:7}  {s['fallback']:8}  {s['no_match']:8}")
print(f"\n✅  CSV written → {out.resolve()}")
