import re
import csv
from pathlib import Path
from collections import defaultdict

logs = {
    "starting_populations_0_3.log": "0.3",
    "starting_populations_0_4.log": "0.4",
    "starting_populations_0_5.log": "0.5",
}

pattern = re.compile(
    r"Number of final molecules for\s+(\S+):\s+(\d+)", re.IGNORECASE
)

table = defaultdict(dict)
order = []   # ← will hold targets in first‐seen order

for fname, column in logs.items():
    path = Path(fname)
    if not path.is_file():
        print(f"⚠️  Skipping {fname}: not found")
        continue

    text = path.read_text()
    for m in pattern.finditer(text):
        target, count = m.groups()
        count = int(count)
        if target not in order:
            order.append(target)
        table[target][column] = count

columns = ["0.3", "0.4", "0.5"]
out_file = Path("Starting_populations.csv")

with out_file.open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["target", *columns])

    for idx, target in enumerate(order, start=1):
        display_name = f"{idx}_{target}"
        row = [display_name] + [table[target].get(col, "") for col in columns]
        writer.writerow(row)

print(f"✅  Wrote {out_file.resolve()}")
