#!/usr/bin/env python3
# A2 step-2 specificity control: m6A D1-HC50 verdict sites x REDIportal overlap.
# m6A sites are likewise A-centred — if ANY A-centred caller output enriched at
# REDIportal positions, the inosine enrichment would be non-specific. Expected:
# m6A overlap ~= background (m6A and A-to-I editing are unrelated chemistries).
# 2026-10-22.

import gzip, re, sys
from collections import defaultdict

ROOT = "/mnt/d/desktop/DRS/DRS_TechShowcase"
DATA = f"{ROOT}/fig3_8mod_spatial/data"
RD = f"{DATA}/rediportal_20261022"
GTF = f"{RD}/ref_grcm38/Mus_musculus.GRCm38.102.gtf.gz"
SL = f"{RD}/sitelevel"
RPT = f"{DATA}/rediportal_20261022/TABLE1_mm10_v3.txt.gz"
M6A = f"{DATA}/site_ctx_nc_day1_NEWPATH.tsv.gz"
REACH = {"exonic", "UTR3", "UTR5", "splicing", "ncRNA_exonic", "ncRNA_splicing"}

# ---- GTF exon chains (same as SiteLevel_02) ----
tid_re = re.compile(r'transcript_id "(ENSMUST\d+)(?:\.\d+)?"')
exons = defaultdict(list); txinfo = {}
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = tid_re.search(f[8])
        if not m: continue
        exons[m.group(1)].append((int(f[3]), int(f[4])))
        txinfo[m.group(1)] = (f[0], f[6])
for t in exons: exons[t].sort()

def lift(t, pos1):
    if t not in exons: return None
    chrom, strand = txinfo[t]; ex = exons[t]; acc = 0
    it = ex if strand == "+" else reversed(ex)
    for s, e in it:
        L = e - s + 1
        if acc + L >= pos1:
            g = s + (pos1 - acc - 1) if strand == "+" else e - (pos1 - acc - 1)
            return ("chr" + ("M" if chrom == "MT" else chrom), g, strand)
        acc += L
    return None

# ---- REDIportal reachable ----
rp = set()
with gzip.open(RPT, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[9] in REACH:
            rp.add((f[0], int(f[1]), f[4]))

# ---- m6A verdict D1-HC50 + paired universe ----
verdict = []; universe = []
with gzip.open(M6A, "rt") as fh:
    head = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[4] == "" or f[5] == "": continue
        nNC, nD1 = int(f[4]), int(f[5])
        if nNC != 3 or nD1 != 3: continue
        if f[6] == "" or f[7] == "": continue
        mS, mD = float(f[6]), float(f[7])
        if not f[0].startswith("ENSMUST"): continue
        universe.append((f[0], int(f[1])))
        if mD >= 0.5:
            verdict.append((f[0], int(f[1])))
print(f"m6A paired ENSMUST universe: {len(universe)}; D1-HC50 verdict: {len(verdict)}")
assert len(verdict) > 0
# canonical all-contig n=13648 — ENSMUST subset expected slightly smaller

def ov(sites):
    n = 0; lifted = 0
    for t, p in sites:
        r = lift(t, p)
        if r is None: continue
        lifted += 1
        if (r[0], r[1], r[2]) in rp: n += 1
    return lifted, n

lU, hU = ov(universe)
lV, hV = ov(verdict)
rU = hU / lU if lU else 0
rV = hV / lV if lV else 0
enr = rV / rU if rU else float("nan")
try:
    from scipy.stats import fisher_exact
    p = fisher_exact([[hV, lV - hV], [hU, lU - hU]])[1]
except Exception:
    p = float("nan")

rep = (
    "m6A cross-modification control x REDIportal (ENSMUST subset)\n"
    + "=" * 60 + "\n"
    f"m6A paired universe: {lU} lifted, {hU} overlap ({100*rU:.4f}%)\n"
    f"m6A D1-HC50 verdict: {lV} lifted, {hV} overlap ({100*rV:.4f}%)\n"
    f"enrichment: {enr:.2f}-fold; Fisher P: {p:.3e}\n"
    "CONTRAST inosine verdict: 5.98% overlap, 114-fold, P=5.2e-13\n"
)
with open(f"{SL}/m6A_control_report.txt", "w") as out:
    out.write(rep)
print(rep)
