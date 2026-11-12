#!/usr/bin/env python3
# LiftChain step-1: discover the company's ENSMUSG<tN> transcript-numbering rule.
# 2026-10-22 — A2 step-2 (site-level REDIportal overlap) prerequisite.
#
# Idea: every diff-site row carries the caller-reported 5-mer motif centred on
# the modified base. For ENSMUSG*****.tN contigs the TRUE transcript sequence
# must reproduce that motif at the site position. Test candidate numbering
# rules against the Ensembl-102 cDNA FASTA (byte-identical to local reference):
#   R1 fasta_order   : Nth transcript of the gene in cDNA FASTA order
#   R2 id_sorted     : Nth transcript sorted by ENSMUST base id
#   R3 len_desc      : Nth transcript by descending sequence length
#   R4 len_asc       : Nth transcript by ascending sequence length
# Position convention: pos 1-based (seq[pos-1]=centre) vs 0-based (seq[pos]=centre).
# The winning rule should reproduce ~100% of motifs (caller-side motif is read
# from the same reference the contigs were built from).

import gzip, re, sys
from collections import defaultdict

DATA = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
FASTA = "/mnt/d/desktop/DRS/DRS_TechShowcase/reference/Mus_musculus.GRCm38.cdna.all.fa.gz"
SITES = f"{DATA}/inosine_diff_site_context_NEWPATH.tsv.gz"
OUT = f"{DATA}/rediportal_20261022/liftchain_rule_discovery.txt"

# ---------- load cDNA FASTA ----------
gene_tx = defaultdict(list)   # ENSMUSG_base -> [(ensmust_base, seq)] in FASTA order
seqs = {}                       # ensmust_base -> seq
hdr_re = re.compile(r"^>(\S+).*gene:(\S+?)\s")
name, gene = None, None
buf = []
def flush():
    if name is None: return
    base_t = name.split(".")[0]
    base_g = gene.split(".")[0]
    s = "".join(buf).upper()
    seqs[base_t] = s
    gene_tx[base_g].append((base_t, s))

with gzip.open(FASTA, "rt") as fh:
    for line in fh:
        if line.startswith(">"):
            flush()
            m = hdr_re.match(line)
            name = m.group(1)
            gene = m.group(2)
            buf = []
        else:
            buf.append(line.strip())
flush()
sys.stderr.write(f"FASTA: {len(seqs)} transcripts, {len(gene_tx)} genes\n")

# ---------- load ENSMUSG.tN sites ----------
sites = []  # (gene_base, tN, pos_int, motif)
with gzip.open(SITES, "rt") as fh:
    head = fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if not f[0].startswith("ENSMUSG"):
            continue
        g, t = f[0].rsplit(".", 1)
        if not t.startswith("t"):
            continue
        sites.append((g, int(t[1:]), int(f[1]), f[7].upper()))
sys.stderr.write(f"ENSMUSG.tN sites: {len(sites)}\n")

# ---------- candidate rules ----------
def rule_fasta_order(gene, n):
    lst = gene_tx.get(gene, [])
    return lst[n-1][0] if 1 <= n <= len(lst) else None

def rule_id_sorted(gene, n):
    lst = sorted(gene_tx.get(gene, []), key=lambda x: x[0])
    return lst[n-1][0] if 1 <= n <= len(lst) else None

def rule_len_desc(gene, n):
    lst = sorted(gene_tx.get(gene, []), key=lambda x: -len(x[1]))
    return lst[n-1][0] if 1 <= n <= len(lst) else None

def rule_len_asc(gene, n):
    lst = sorted(gene_tx.get(gene, []), key=lambda x: len(x[1]))
    return lst[n-1][0] if 1 <= n <= len(lst) else None

rules = {"R1_fasta_order": rule_fasta_order, "R2_id_sorted": rule_id_sorted,
         "R3_len_desc": rule_len_desc, "R4_len_asc": rule_len_asc}

# ---------- score ----------
res = {r: {"match1": 0, "match0": 0, "centreA1": 0, "centreA0": 0, "tot": 0, "miss_gene": 0, "miss_n": 0}
       for r in rules}
example_map = {}  # (gene,tN) -> chosen ENSMUST under winning rule, filled later
for g, n, pos, motif in sites:
    if g not in gene_tx:
        for r in rules: res[r]["miss_gene"] += 1
        continue
    for rname, rfun in rules.items():
        t = rfun(g, n)
        R = res[rname]
        if t is None:
            R["miss_n"] += 1
            continue
        R["tot"] += 1
        s = seqs[t]
        # 1-based: centre at index pos-1, 5-mer = s[pos-3:pos+2]
        m1 = s[pos-3:pos+2] if 2 <= pos-1 < len(s)-1 else ""
        m0 = s[pos-2:pos+3] if 1 <= pos < len(s)-2 else ""
        if m1 == motif: R["match1"] += 1
        if m0 == motif: R["match0"] += 1
        if len(m1) == 5 and m1[2] == "A": R["centreA1"] += 1
        if len(m0) == 5 and m0[2] == "A": R["centreA0"] += 1

with open(OUT, "w") as out:
    out.write(f"# tN numbering-rule discovery — {len(sites)} ENSMUSG.tN inosine diff sites\n")
    out.write("# match1/match0 = exact 5-mer motif reproduction (1-based / 0-based pos)\n")
    out.write("# centreA = centre base is A (caller invariant for inosine)\n\n")
    out.write("rule\ttestable\tmotif_match_1based\tmotif_match_0based\tcentreA_1based\tcentreA_0based\tgene_missing\tn_out_of_range\n")
    for rname in rules:
        R = res[rname]
        out.write(f"{rname}\t{R['tot']}\t{R['match1']}\t{R['match0']}\t{R['centreA1']}\t{R['centreA0']}\t{R['miss_gene']}\t{R['miss_n']}\n")
    out.write("\n# gene universe sizes\n")
    out.write(f"# genes in FASTA: {len(gene_tx)}; transcripts: {len(seqs)}\n")
print(open(OUT).read())
