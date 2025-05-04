import re
import csv
from pathlib import Path
from collections import defaultdict

logs = {
    "log_pops_0_3.txt": "0.3",
    "log_pops_0_4.txt": "0.4",
    "log_pops_0_5.txt": "0.5",      
}

pattern = re.compile(
    r"Number of final molecules for\s+(\S+):\s+(\d+)", re.IGNORECASE
)

table = defaultdict(dict)      

for fname, column in logs.items():
    path = Path(fname)
    if not path.is_file():
        print(f"⚠️  Skipping {fname}: not found")
        continue

    with path.open() as f:
        for m in pattern.finditer(f.read()):
            target, count = m.groups()
            table[target][column] = int(count)


columns = ["0.3", "0.4", "0.5"] 
out_file = Path("Starting_populations.csv")

with out_file.open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["target", *columns])

    for target in sorted(table):
        row = [target] + [table[target].get(col, "") for col in columns]
        writer.writerow(row)

print(f"✅  Wrote {out_file.resolve()}")
