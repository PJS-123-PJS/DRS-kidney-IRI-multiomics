#!/usr/bin/env python3
# LiftChain step-2: test whether ENSMUSG*****.tN contigs follow the MGP GTF's
# per-gene transcript order (mapped to projection-parent ENSMUST sequence).
# 2026-10-22 — A2 step-2 prerequisite, runs fully on local D: files.
#
# Outputs a compact parent-mapping TSV for the lift chain + rule test stats.

import gzip, re, sys
from collections import defaultdict

GTF = "/mnt/d/desktop/DRS/ATAC-03.align_location/Mus_musculus.gtf"
FASTA = "/mnt/d/desktop/DRS/DRS_TechShowcase/reference/Mus_musculus.GRCm38.cdna.all.fa.gz"
SITES = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/inosine_diff_site_context_NEWPATH.tsv.gz"
OUTDIR = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/rediportal_20261022"

# ---------- 1. parse MGP GTF transcript lines ----------
attr_re = {
    "gid":  re.compile(r'gene_id "(MGP_C57BL6NJ_G\d+)"'),
    "tid":  re.compile(r'transcript_id "(MGP_C57BL6NJ_T\d+)"'),
    "pg":   re.compile(r'projection_parent_gene "(ENSMUSG\d+)\.\d+"'),
    "pt":   re.compile(r'projection_parent_transcript "(ENSMUST\d+)\.\d+"'),
}
parent_tx = defaultdict(list)   # ENSMUSG_base -> [(mgp_t, parent_ensmust_base)] in GTF order
mgp_info = {}                   # mgp_t -> (chrom, tstart, tend, strand, parent_g, parent_t)
n_lines = 0
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript": continue
        n_lines += 1
        a = f[8]
        m_gid = attr_re["gid"].search(a)
        m_tid = attr_re["tid"].search(a)
        m_pg  = attr_re["pg"].search(a)
        m_pt  = attr_re["pt"].search(a)
        if not (m_gid and m_tid): continue
        mgp_t = m_tid.group(1)
        pg = m_pg.group(1) if m_pg else ""
        pt = m_pt.group(1) if m_pt else ""
        mgp_info[mgp_t] = (f[0], int(f[3]), int(f[4]), f[6], pg, pt)
        if pg:
            parent_tx[pg].append((mgp_t, pt))
sys.stderr.write(f"GTF transcripts parsed: {n_lines}; parent genes: {len(parent_tx)}\n")

# ---------- 2. load cDNA sequences (release-102, parent ENSMUST) ----------
seqs = {}
name = None; buf = []
with gzip.open(FASTA, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            if name: seqs[name.split(".")[0]] = "".join(buf).upper()
            name = line.split()[0][1:]; buf = []
        else: buf.append(line.strip())
if name: seqs[name.split(".")[0]] = "".join(buf).upper()
sys.stderr.write(f"cDNA transcripts: {len(seqs)}\n")

# ---------- 3. load ENSMUSG.tN sites ----------
sites = []
with gzip.open(SITES, "rt") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if not f[0].startswith("ENSMUSG"): continue
        g, t = f[0].rsplit(".", 1)
        sites.append((g, int(t[1:]), int(f[1]), f[7].upper()))

# ---------- 4. test rules ----------
def r_gtf_order(g, n):
    lst = parent_tx.get(g, [])
    return lst[n-1] if 1 <= n <= len(lst) else None
def r_mgp_sorted(g, n):
    lst = sorted(parent_tx.get(g, []), key=lambda x: x[0])
    return lst[n-1] if 1 <= n <= len(lst) else None
def r_parent_sorted(g, n):
    lst = sorted(parent_tx.get(g, []), key=lambda x: x[1])
    return lst[n-1] if 1 <= n <= len(lst) else None
rules = {"M1_gtf_order": r_gtf_order, "M2_mgp_sorted": r_mgp_sorted, "M3_parent_sorted": r_parent_sorted}
res = {r: dict(match=0, tot=0, no_seq=0, miss_g=0, miss_n=0) for r in rules}
for g, n, pos, motif in sites:
    if g not in parent_tx:
        for r in rules: res[r]["miss_g"] += 1
        continue
    for rn, rf in rules.items():
        R = res[rn]
        pick = rf(g, n)
        if pick is None: R["miss_n"] += 1; continue
        pt = pick[1]
        if pt not in seqs: R["no_seq"] += 1; continue
        s = seqs[pt]
        R["tot"] += 1
        m1 = s[pos-3:pos+2] if 2 <= pos-1 < len(s)-1 else ""
        if m1 == motif: R["match"] += 1

with open(f"{OUTDIR}/liftchain_rule_discovery_MGP.txt", "w") as out:
    out.write(f"# MGP-order rule test — {len(sites)} ENSMUSG.tN sites, 1-based pos, parent ENSMUST 6J sequence\n")
    out.write("rule\ttestable\tmotif_match\tmatch_pct\tgene_missing\tn_out_of_range\tparent_seq_missing\n")
    for rn in rules:
        R = res[rn]
        pct = 100.0*R["match"]/R["tot"] if R["tot"] else 0
        out.write(f"{rn}\t{R['tot']}\t{R['match']}\t{pct:.2f}\t{R['miss_g']}\t{R['miss_n']}\t{R['no_seq']}\n")
print(open(f"{OUTDIR}/liftchain_rule_discovery_MGP.txt").read())

# ---------- 5. dump parent mapping table for the lift chain ----------
with open(f"{OUTDIR}/mgp_parent_transcripts.tsv", "w") as out:
    out.write("parent_gene\torder_in_gtf\tmgp_transcript\tparent_transcript\tchrom\ttstart\ttend\tstrand\n")
    for g, lst in parent_tx.items():
        for i, (mt, pt) in enumerate(lst, 1):
            c, s_, e_, sd, _, _ = mgp_info[mt]
            out.write(f"{g}\t{i}\t{mt}\t{pt}\t{c}\t{s_}\t{e_}\t{sd}\n")
sys.stderr.write("mapping table written\n")
