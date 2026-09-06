"""Extract FinnGen DR outcome data for the 25 drug-target cis-IVs (streaming filter)."""
import gzip, csv

RSID_FILE = "E:/workbuddy/2026-09-04-18-19-16/phase1_out/iv_druggable.txt"
GZ = "E:/workbuddy/2026-09-04-18-19-16/phase1_out/finngen_R13_DM_RETINOPATHY.gz"
OUT = "E:/workbuddy/2026-09-04-18-19-16/phase1_out/finngen_drug_subset.tsv"

want = set()
for line in open(RSID_FILE):
    want.add(line.strip())
print("target rsids:", len(want))

hit = 0
with gzip.open(GZ, "rt", encoding="utf-8", errors="replace") as f, \
     open(OUT, "w", newline="") as out:
    w = csv.writer(out, delimiter="\t", lineterminator="\n")
    hdr = f.readline().rstrip("\n")
    cols = hdr.lstrip("#").split("\t")
    w.writerow(cols)
    ix = {c: i for i, c in enumerate(cols)}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) <= 5:
            continue
        rs = parts[ix["rsids"]]
        if rs in want:
            w.writerow(parts)
            hit += 1
            if hit == len(want):
                break
print(f"matched {hit}/{len(want)} rsids -> {OUT}")
