"""
Parse CAE matching log files (with INFO-level lines) into a summary CSV.
Usage: python parse_cae_matching_log.py cae_analysis_0_5.log cae_summary_0_5.csv
"""

import sys, re, csv
from pathlib import Path

# --- 1. CLI & file check ---
if len(sys.argv) != 3:
    sys.exit("usage: python parse_matching_log.py <logfile> <output.csv>")

log_path = Path(sys.argv[1])
if not log_path.is_file():
    sys.exit(f"file not found: {log_path}")

out = Path(sys.argv[2])

# --- 2. Regexes for new log style ---
re_target    = re.compile(r"INFO:\s*🔍 Target:\s+(\S+)")
re_matched   = re.compile(r"INFO:\s*📊 Matched\s+(\d+)/\d+\s+CAEs")
re_distorted = re.compile(r"INFO:\s*🟡 (\d+)\s+matched only as Distorted Matches")
re_nomatch   = re.compile(r"INFO:\s*🟥 (\d+)\s+had no population CAE")

summary, current = {}, None
stats = {"matched": 0, "distorted": 0, "no_match": 0}

with log_path.open() as fh:
    for line in fh:
        if m := re_target.search(line):
            if current is not None:
                summary[current] = stats
                stats = {"matched": 0, "distorted": 0, "no_match": 0}
            current = m.group(1)
        elif m := re_matched.search(line):
            stats["matched"] = int(m.group(1))
        elif m := re_distorted.search(line):
            stats["distorted"] = int(m.group(1))
        elif m := re_nomatch.search(line):
            stats["no_match"] = int(m.group(1))

if current is not None:  # Flush last
    summary[current] = stats

# --- 3. Write CSV ---
with out.open("w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["target", "matched", "distorted", "no_match"])
    for tgt, s in summary.items():
        w.writerow([tgt, s["matched"], s["distorted"], s["no_match"]])

# --- 4. Console preview ---
print(f"📄  Results for {log_path.name}")
print("target          matched  distorted  no_match")
print("──────────────  ───────  ─────────  ────────")
for tgt, s in summary.items():
    print(f"{tgt:14}  {s['matched']:7}  {s['distorted']:9}  {s['no_match']:8}")
print(f"\n✅  CSV written → {out.resolve()}")
