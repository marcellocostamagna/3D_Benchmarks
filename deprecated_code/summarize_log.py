import re
from collections import defaultdict

def summarize_log(log_text):
    targets = []
    current_target = None
    total_fragments = 0
    matched_fragments = 0
    fragment_results = []

    for line in log_text.splitlines():
        # New target
        if match := re.search(r"🔍 Target: (\w+)", line):
            if current_target:
                targets.append({
                    "name": current_target,
                    "total": total_fragments,
                    "matched": matched_fragments,
                    "fragments": fragment_results
                })
            current_target = match.group(1)
            total_fragments = 0
            matched_fragments = 0
            fragment_results = []

        # Fragment count
        if match := re.search(r"✅ Unique target fragments: (\d+)", line):
            total_fragments = int(match.group(1))

        # Matched fragment
        if match := re.search(r"🔎 Matched (\d+)/(\d+)", line):
            matched, total = int(match.group(1)), int(match.group(2))
            fragment_results.append((matched, total))

        # Summary matched count
        if match := re.search(r"📊 Matched (\d+)/(\d+)", line):
            matched_fragments = int(match.group(1))

    # Append last target
    if current_target:
        targets.append({
            "name": current_target,
            "total": total_fragments,
            "matched": matched_fragments,
            "fragments": fragment_results
        })

    # Print summary
    for idx, t in enumerate(targets, 1):
        mismatched = []
        sum_matches = sum(m for m, _ in t["fragments"])
        sum_expected = sum(n for _, n in t["fragments"])
        if sum_matches != t["matched"] or t["total"] != len(t["fragments"]):
            # Find which fragments had mismatches
            for i, (m, n) in enumerate(t["fragments"], 1):
                if m != n:
                    mismatched.append(f"frag {i} ({m}/{n})")
            print(f"{idx}_{t['name']}: {t['total']} fragments, MISMATCH – {t['matched']}/{t['total']} claimed matched, but mismatch in: {', '.join(mismatched) or 'unknown'}")
        else:
            print(f"{idx}_{t['name']}: {t['total']} fragments, ALL Matched!")

# Example usage:
with open("log_matching_0_3.txt", "r") as f:
    log_data = f.read()
summarize_log(log_data)
