#!/usr/bin/env python3
# A->G mismatch hard-signature: build target site lists (inosine verdict/diffUp
# + m6A DRACH+ SH-HC70 control) as 0-based BED keys for the modkit extraction.
# 2026-10-22.

import gzip, os

ROOT = "/mnt/d/desktop/DRS/DRS_TechShowcase"
DATA = f"{ROOT}/fig3_8mod_spatial/data"
SL = f"{DATA}/rediportal_20261022/sitelevel"
OUTD = f"{SL}/mismatch"
os.makedirs(OUTD, exist_ok=True)

rows = []

# inosine verdict (141) — pos is 1-based -> start0 = pos-1
with open(f"{SL}/verdict_D1HC50_sites.tsv") as fh:
    fh.readline()
    for line in fh:
        t, p = line.rstrip("\n").split("\t")[:2]
        rows.append((t, int(p) - 1, "ino_verdict"))

# inosine diff Up (11,233)
with open(f"{SL}/diff_up_sites.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        rows.append((f[0], int(f[1]) - 1, "ino_diffup"))

# m6A DRACH+ SH-HC70 control (canonical n=5,174): paired 3/3, is_drach, mean_frac_Sham>=0.7
n_m6a = 0
with gzip.open(f"{DATA}/site_ctx_nc_day1_NEWPATH.tsv.gz", "rt") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) < 8 or f[4] == "" or f[5] == "" or f[6] == "" or f[7] == "":
            continue
        if int(f[4]) != 3 or int(f[5]) != 3:
            continue
        if f[2] != "TRUE":
            continue
        if float(f[6]) >= 0.7:
            rows.append((f[0], int(f[1]) - 1, "m6A_ctrl"))
            n_m6a += 1
print(f"m6A control sites: {n_m6a} (canonical 5174)")
assert n_m6a == 5174, f"m6A control count mismatch: {n_m6a}"

with open(f"{OUTD}/target_sites.bedkeys.tsv", "w") as out:
    for t, s0, cls in rows:
        out.write(f"{t}\t{s0}\t{cls}\n")
print(f"total target keys: {len(rows)}")
