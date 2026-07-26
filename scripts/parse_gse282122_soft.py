import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path

source = Path("data/metadata/GSE282122_family.soft.gz")
output = Path("results/feasibility/GSE282122_sample_metadata.csv")

records = []
record = None
with gzip.open(source, "rt", encoding="utf-8", errors="replace") as handle:
    for raw in handle:
        line = raw.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if record:
                records.append(record)
            record = {"accession": line.split(" = ", 1)[1]}
        elif record is not None and line.startswith("!Sample_title = "):
            record["title"] = line.split(" = ", 1)[1]
        elif record is not None and line.startswith("!Sample_characteristics_ch1 = "):
            value = line.split(" = ", 1)[1]
            if ": " in value:
                key, val = value.split(": ", 1)
                record[key.lower()] = val
if record:
    records.append(record)

fields = sorted({key for row in records for key in row})
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w", newline="", encoding="utf-8-sig") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields)
    writer.writeheader()
    writer.writerows(records)

print(f"samples={len(records)}")
for field in ["disease", "treatment", "inflammation", "site", "batch", "match", "librarytype"]:
    print(field, dict(Counter(row.get(field, "MISSING") for row in records)))

patients = defaultdict(list)
for row in records:
    patients[row.get("patient", "MISSING")].append(row)
print(f"patients={len(patients)}")
for disease in sorted({r.get('disease', 'MISSING') for r in records}):
    ids = {r.get("patient") for r in records if r.get("disease") == disease}
    paired = {
        pid for pid in ids
        if {r.get("treatment") for r in patients[pid]} >= {"Pre", "Post"}
    }
    print(f"{disease}: patients={len(ids)}, pre_post_patients={len(paired)}")

keys = sorted({key for row in records for key in row})
outcome_like = [key for key in keys if any(term in key for term in ("remission", "response", "outcome", "mayo", "harvey", "cda"))]
print("outcome_fields", outcome_like)
