#!/usr/bin/env python3
# A2 step-2 core: lift inosine sites (ENSMUST transcript-relative, 1-based) to
# mm10 genomic coordinates via Ensembl-102 GRCm38 GTF exon chains, then exact
# position overlap vs REDIportal TABLE1_mm10_v3 (transcriptome-reachable subset).
# 2026-10-22.
#
# Inputs
#   ref_grcm38/Mus_musculus.GRCm38.102.gtf.gz   (exon chains)
#   reference/Mus_musculus.GRCm38.cdna.all.fa.gz (A-position null universe)
#   sitelevel/{verdict_D1HC50_sites,diff_up_sites}.tsv, detected_universe.tsv.gz
#   TABLE1_mm10_v3.txt.gz (REDIportal)
# Outputs (sitelevel/)
#   lifted_{verdict,diffup,universe}.bed.gz     — chrom,start0,end,strand,trans,pos
#   overlap_report.txt                          — overlap counts + OR + Fisher
#   overlap_hits.tsv                            — per-site hit table
#
# Lift conventions: site pos is 1-based on spliced transcript (validated by
# 99.4% motif reproduction). Genomic coordinate = exon-walk, strand-aware.
# REDIportal chrom = chrN; GTF chrom = N. MT -> chrM mapping applied.

import gzip, re, sys
from collections import defaultdict

ROOT = "/mnt/d/desktop/DRS/DRS_TechShowcase"
DATA = f"{ROOT}/fig3_8mod_spatial/data"
RD = f"{DATA}/rediportal_20261022"
GTF = f"{RD}/ref_grcm38/Mus_musculus.GRCm38.102.gtf.gz"
FASTA = f"{ROOT}/reference/Mus_musculus.GRCm38.cdna.all.fa.gz"
SL = f"{RD}/sitelevel"
RPT = f"{DATA}/rediportal_20261022/TABLE1_mm10_v3.txt.gz"

REACH = {"exonic", "UTR3", "UTR5", "splicing", "ncRNA_exonic", "ncRNA_splicing"}

# ---------- 1. GTF: exon chains per transcript ----------
sys.stderr.write("parse GTF exons...\n")
tid_re = re.compile(r'transcript_id "(ENSMUST\d+)(?:\.\d+)?"')
exons = defaultdict(list)   # ensmust_base -> [(gstart1, gend1)]  (1-based inclusive, GTF native)
txinfo = {}                 # ensmust_base -> (chrom, strand)
with gzip.open(GTF, "rt") as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = tid_re.search(f[8])
        if not m: continue
        t = m.group(1)
        exons[t].append((int(f[3]), int(f[4])))
        txinfo[t] = (f[0], f[6])
for t in exons:
    exons[t].sort()   # genomic order; sense order handled at lift time
sys.stderr.write(f"transcripts with exons: {len(exons)}\n")

def lift(t, pos1):
    """transcript-relative 1-based pos -> (chrom, gpos1, strand) or None"""
    if t not in exons: return None
    chrom, strand = txinfo[t]
    ex = exons[t]
    acc = 0
    if strand == "+":
        for s, e in ex:
            L = e - s + 1
            if acc + L >= pos1:
                return (chrom, s + (pos1 - acc - 1), strand)
            acc += L
    else:
        for s, e in reversed(ex):
            L = e - s + 1
            if acc + L >= pos1:
                return (chrom, e - (pos1 - acc - 1), strand)
            acc += L
    return None

# ---------- 2. load site sets ----------
def load_sites(path, filt):
    out = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if filt(f[0]):
                out.append((f[0], int(f[1])))
    return out

ensmust_only = lambda t: t.startswith("ENSMUST")
verdict = load_sites(f"{SL}/verdict_D1HC50_sites.tsv", ensmust_only)
diffup  = load_sites(f"{SL}/diff_up_sites.tsv", ensmust_only)
universe = []
with gzip.open(f"{SL}/detected_universe.tsv.gz", "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if ensmust_only(f[0]):
            universe.append((f[0], int(f[1])))
sys.stderr.write(f"sets: verdict {len(verdict)}, diffup {len(diffup)}, universe {len(universe)}\n")

# ---------- 3. lift ----------
def lift_all(sites, outname):
    bed = []
    fail = 0
    for t, p in sites:
        r = lift(t, p)
        if r is None: fail += 1; continue
        chrom, gpos, strand = r
        c = "chr" + ("M" if chrom == "MT" else chrom)
        bed.append((c, gpos - 1, gpos, strand, t, p))
    with gzip.open(f"{SL}/{outname}", "wt") as out:
        for b in sorted(bed):
            out.write("\t".join(map(str, b)) + "\n")
    sys.stderr.write(f"{outname}: lifted {len(bed)}, failed {fail}\n")
    return bed

bed_v = lift_all(verdict, "lifted_verdict.bed.gz")
bed_d = lift_all(diffup, "lifted_diffup.bed.gz")
bed_u = lift_all(universe, "lifted_universe.bed.gz")

# ---------- 4. REDIportal reachable subset ----------
rp = {}          # (chrom, pos1, strand) -> gene
rp_anystrand = set()
with gzip.open(RPT, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[9] not in REACH: continue
        key = (f[0], int(f[1]), f[4])
        rp[key] = f[10]
        rp_anystrand.add((f[0], int(f[1])))
sys.stderr.write(f"REDIportal reachable sites: {len(rp)} (any-strand {len(rp_anystrand)})\n")

# ---------- 5. overlap ----------
def overlap(bed):
    hits, hits_any = [], 0
    for c, s0, e, strand, t, p in bed:
        key = (c, e, strand)
        if key in rp:
            hits.append((c, e, strand, t, p, rp[key]))
        elif (c, e) in rp_anystrand:
            hits_any += 1
    return hits, hits_any

hv, hav = overlap(bed_v)
hd, had = overlap(bed_d)
hu, hau = overlap(bed_u)

# ---------- 6. null: fraction of A positions on detected transcripts in REDIportal ----------
sys.stderr.write("compute A-position null...\n")
det_tx = set(t for t, _ in universe)
hdr_re = re.compile(r"^>(\S+)")
A_tot = 0
with gzip.open(FASTA, "rt") as fh:
    name = None; buf = []
    def flush2():
        global A_tot
        if name is None: return
        base = name.split(".")[0]
        if base in det_tx:
            A_tot += "".join(buf).upper().count("A")
    for line in fh:
        if line.startswith(">"):
            flush2()
            name = hdr_re.match(line).group(1); buf = []
        else:
            buf.append(line.strip())
    flush2()
sys.stderr.write(f"detected transcripts: {len(det_tx)}, A positions: {A_tot}\n")

# ---------- 7. report ----------
def fisher(a, b, c, d):
    from math import lgamma, exp
    # hypergeometric two-sided via scipy if available else one-sided approx
    try:
        from scipy.stats import fisher_exact
        return fisher_exact([[a, b], [c, d]])[1]
    except Exception:
        return float("nan")

nU, nV, nD = len(bed_u), len(bed_v), len(bed_d)
hU, hV, hD = len(hu), len(hv), len(hd)
base_rate = hU / nU if nU else 0
with open(f"{SL}/overlap_report.txt", "w") as out:
    out.write("A2 step-2 site-level REDIportal overlap (ENSMUST subset, exact mm10 lift)\n")
    out.write("=" * 76 + "\n")
    out.write(f"REDIportal transcriptome-reachable sites: {len(rp)}\n")
    out.write(f"Universe (detected ENSMUST sites): {nU}, overlap: {hU} ({100*hU/nU:.4f}%)\n")
    out.write(f"Verdict D1-HC50 (ENSMUST): {nV}, overlap: {hV} ({100*hV/max(nV,1):.2f}%)\n")
    out.write(f"Diff-Up (ENSMUST): {nD}, overlap: {hD} ({100*hD/max(nD,1):.2f}%)\n")
    out.write(f"Opposite-strand/basematch-only extra hits: verdict {hav}, diffup {had}, universe {hau}\n")
    if nV and hV:
        out.write(f"Verdict enrichment vs universe: {hV/nV/base_rate:.1f}-fold\n")
    if nD and hD:
        out.write(f"DiffUp enrichment vs universe: {hD/nD/base_rate:.1f}-fold\n")
    pv = fisher(hV, nV - hV, hU, nU - hU)
    pd = fisher(hD, nD - hD, hU, nU - hU)
    out.write(f"Fisher P verdict vs universe: {pv:.3e}\n")
    out.write(f"Fisher P diffup vs universe: {pd:.3e}\n")
    out.write(f"A positions on detected transcripts (null denom): {A_tot}\n")

with open(f"{SL}/overlap_hits.tsv", "w") as out:
    out.write("set\tchrom\tgpos\tstrand\ttrans\ttrans_pos\trp_gene\n")
    for h in hv:
        out.write("verdict\t" + "\t".join(map(str, h)) + "\n")
    for h in hd:
        out.write("diffup\t" + "\t".join(map(str, h)) + "\n")

print(open(f"{SL}/overlap_report.txt").read())
