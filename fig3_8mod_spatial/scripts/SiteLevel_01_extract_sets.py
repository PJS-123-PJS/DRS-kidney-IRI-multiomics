#!/usr/bin/env python3
# A2 step-2 builder-1: reconstruct exact inosine verdict/diff site sets.
# 2026-10-22 — mirrors ModVerdict_Paired_Stats.R logic on local ctx chunks.
#
# Outputs (rediportal_20261022/sitelevel/):
#   verdict_D1HC50_sites.tsv   — the exact 141 verdict sites (trans,pos,...)
#   diff_up_sites.tsv          — Up rows of the differential-site table
#   detected_universe.tsv.gz   — union of detected (trans,pos) over 6 samples
# Hard assertion: verdict set must be n=141 (canonical modverdict value).

import gzip, os
from collections import defaultdict

DATA = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
OUTD = f"{DATA}/rediportal_20261022/sitelevel"
os.makedirs(OUTD, exist_ok=True)

NC = ["NCs1", "NCs2", "NCs3"]
D1 = ["Day1s1", "Day1s2", "Day1s3"]

# ---------- aggregate ctx chunks (mirror of site_ctx aggregation) ----------
agg = defaultdict(lambda: [0, 0, 0.0, 0.0])  # (trans,pos) -> nNC, nD1, sumNC, sumD1
for s in NC + D1:
    grp = 0 if s in NC else 1
    with gzip.open(f"{DATA}/site_ctx_chunks_inosine/ctx_{s}.tsv.gz", "rt") as fh:
        head = fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]))
            a = agg[key]
            a[grp] += 1
            a[grp + 2] += float(f[2])
print(f"union detected sites: {len(agg)}")

paired = {}
for k, (nnc, nd1, snc, sd1) in agg.items():
    if nnc == 3 and nd1 == 3:
        paired[k] = (snc / 3, sd1 / 3)
print(f"fully paired sites: {len(paired)}")

verdict = {k: v for k, v in paired.items() if v[1] >= 0.5}
hc70 = {k: v for k, v in paired.items() if v[1] >= 0.7}
print(f"D1-HC50 verdict sites: {len(verdict)}  (canonical n=141)")
print(f"D1-HC70 sites: {len(hc70)}  (canonical n=23)")
assert len(verdict) == 141, f"VERDICT COUNT MISMATCH: {len(verdict)} != 141"
assert len(hc70) == 23, f"HC70 COUNT MISMATCH: {len(hc70)} != 23"

with open(f"{OUTD}/verdict_D1HC50_sites.tsv", "w") as out:
    out.write("trans\tpos\tmean_frac_Sham\tmean_frac_Day1\n")
    for (t, p), (mnc, md1) in sorted(verdict.items()):
        out.write(f"{t}\t{p}\t{mnc:.6f}\t{md1:.6f}\n")

# ---------- diff Up sites ----------
n_up = 0
with gzip.open(f"{DATA}/inosine_diff_site_context_NEWPATH.tsv.gz", "rt") as fh, \
     open(f"{OUTD}/diff_up_sites.tsv", "w") as out:
    out.write("trans\tpos\tname\ta_pct\tb_pct\teffect_size\tmotif\tgene_struc\n")
    head = fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[3] != "Up":
            continue
        n_up += 1
        out.write(f"{f[0]}\t{f[1]}\t{f[2]}\t{f[4]}\t{f[5]}\t{f[6]}\t{f[7]}\t{f[9]}\n")
print(f"diff Up sites: {n_up}")

# ---------- detected universe ----------
with gzip.open(f"{OUTD}/detected_universe.tsv.gz", "wt") as out:
    out.write("trans\tpos\tn_det\tmean_frac_all\n")
    for (t, p), (nnc, nd1, snc, sd1) in sorted(agg.items()):
        out.write(f"{t}\t{p}\t{nnc+nd1}\t{(snc+sd1)/(nnc+nd1):.6f}\n")
print("DONE")
