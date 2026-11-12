#!/usr/bin/env python3
# A->G mismatch hard-signature analysis (modkit Ndiff-based).
# 2026-10-22 — A2 step-2 decisive arm.
#
# Logic: inosine is mis-basecalled as G at some rate; such reads land in
# Ndiff (base differs from reference A). True inosine sites should show
# elevated mismatch share, scaling with mod fraction, and rising NC->Day1
# at sites that gain inosine. m6A HC sites (A-centred, different chemistry)
# are the specificity control.
#
# mismatch share = Ndiff / (Nvalid + Ndiff)   [confident basecalls]

import os, sys
from collections import defaultdict

ROOT = "/mnt/d/desktop/DRS/DRS_TechShowcase"
SL = f"{ROOT}/fig3_8mod_spatial/data/rediportal_20261022/sitelevel"
MM = f"{SL}/mismatch"
NC = ["NCs1", "NCs2", "NCs3"]
D1 = ["Day1s1", "Day1s2", "Day1s3"]

# ---- membership sets (re-derived, overwrite-proof) ----
verdict = set()
with open(f"{SL}/verdict_D1HC50_sites.tsv") as fh:
    fh.readline()
    for line in fh:
        t, p = line.rstrip("\n").split("\t")[:2]
        verdict.add((t, int(p) - 1))
diffup = set()
with open(f"{SL}/diff_up_sites.tsv") as fh:
    fh.readline()
    for line in fh:
        f = line.rstrip("\n").split("\t")
        diffup.add((f[0], int(f[1]) - 1))
m6actrl = set()
with open(f"{MM}/target_sites.bedkeys.tsv") as fh:
    for line in fh:
        t, s0, cls = line.rstrip("\n").split("\t")
        if cls == "m6A_ctrl":
            m6actrl.add((t, int(s0)))

# ---- load extracts ----
# per-sample dicts: key -> (Nvalid, Nmod, Ncan, Nfail, Ndiff)
def load(path):
    d = {}
    bg = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 19:
                continue
            nv, nmod, ncan = int(f[13]), int(f[15]), int(f[16])
            nfail, ndiff = int(f[19-3]), int(f[19-2])
            # cols: 0=class,1..18=bed (bed col10=f[10+... let me be explicit:
            # f[0]=class; f[1]=chrom f[2]=start f[3]=end f[4]=name f[5]=score
            # f[6]=strand f[7]=thick f[8]=thick f[9]=color f[10]=Nvalid
            # f[11]=frac f[12]=Nmod f[13]=Ncan f[14]=Nother f[15]=Ndelete
            # f[16]=Nfail f[17]=Ndiff f[18]=Nnocall
            nv, nmod, ncan = int(f[10]), int(f[12]), int(f[13])
            nfail, ndiff = int(f[16]), int(f[17])
            rec = (nv, nmod, ncan, nfail, ndiff)
            if f[0] == "bg":
                bg.append(rec)
            else:
                d[(f[1], int(f[2]))] = rec
    return d, bg

data, bgs = {}, {}
for mod in ["inosine", "m6a"]:
    for s in NC + D1:
        d, b = load(f"{MM}/extract_{mod}_{s}.tsv")
        data[(mod, s)] = d
        bgs[(mod, s)] = b
        print(f"{mod}/{s}: targets={len(d)} bg={len(b)}", file=sys.stderr)

from scipy.stats import mannwhitneyu, fisher_exact, wilcoxon

def share(rec):
    nv, nmod, ncan, nfail, ndiff = rec
    den = nv + ndiff
    return ndiff / den if den > 0 else None

def pool(recs):
    nd = sum(r[4] for r in recs)
    tot = sum(r[0] + r[4] for r in recs)
    return nd, tot

report = []
def w(s=""):
    report.append(s)
    print(s)

# ============ A. verdict sites: NC vs Day1 (per-sample + site-paired) ============
w("== A. verdict sites mismatch share (inosine) ==")
for grp, samples in [("NC", NC), ("Day1", D1)]:
    shares = []
    for s in samples:
        for k in verdict:
            r = data[("inosine", s)].get(k)
            if r:
                sh = share(r)
                if sh is not None:
                    shares.append(sh)
    w(f"{grp}: n_site_sample={len(shares)} mean_share={100*sum(shares)/len(shares):.3f}% "
      f"sites_with_any_mismatch={100*sum(1 for x in shares if x>0)/len(shares):.1f}%")
# pooled counts
for grp, samples in [("NC", NC), ("Day1", D1)]:
    recs = [data[("inosine", s)][k] for s in samples for k in verdict if k in data[("inosine", s)]]
    nd, tot = pool(recs)
    w(f"{grp} pooled: Ndiff={nd} / {tot} = {100*nd/tot:.3f}%")

# ============ B. backgrounds ============
w("\n== B. background A-position mismatch share ==")
for mod in ["inosine", "m6a"]:
    for grp, samples in [("NC", NC), ("Day1", D1)]:
        recs = [r for s in samples for r in bgs[(mod, s)]]
        nd, tot = pool(recs)
        w(f"{mod} bg {grp}: n={len(recs)} pooled Ndiff share={100*nd/tot:.3f}%")

# ============ C. m6A HC control sites ============
w("\n== C. m6A DRACH+ HC70 control sites ==")
for grp, samples in [("NC", NC), ("Day1", D1)]:
    recs = [data[("m6a", s)][k] for s in samples for k in m6actrl if k in data[("m6a", s)]]
    nd, tot = pool(recs)
    shs = [share(r) for r in recs if share(r) is not None]
    w(f"m6A ctrl {grp}: n={len(recs)} pooled share={100*nd/tot:.3f}% mean={100*sum(shs)/len(shs):.3f}%")

# ============ D. dose-response: mismatch share by mod-fraction bin ============
w("\n== D. dose-response (pooled 6 samples, inosine targets=diffup+verdict, m6A=ctrl) ==")
def dose(mod, keys):
    bins = defaultdict(list)
    for s in NC + D1:
        d = data[(mod, s)]
        for k in keys:
            r = d.get(k)
            if not r or r[0] < 10:
                continue
            frac = r[1] / r[0]           # Nmod/Nvalid
            sh = share(r)
            if sh is None: continue
            b = min(int(frac * 5), 4)    # [0,.2)..[.8,1]
            bins[b].append(sh)
    labels = ["0-20%", "20-40%", "40-60%", "60-80%", "80-100%"]
    for b in range(5):
        v = bins.get(b, [])
        if v:
            w(f"  frac {labels[b]}: n={len(v)} mean mismatch={100*sum(v)/len(v):.3f}%")
w("inosine sites:")
dose("inosine", diffup | verdict)
w("m6A ctrl sites:")
dose("m6a", m6actrl)

# ============ E. stats ============
w("\n== E. stats ==")
v_rec_d1 = [data[("inosine", s)][k] for s in D1 for k in verdict if k in data[("inosine", s)]]
v_rec_nc = [data[("inosine", s)][k] for s in NC for k in verdict if k in data[("inosine", s)]]
bg_rec = [r for s in NC + D1 for r in bgs[("inosine", s)]]
m6_rec = [data[("m6a", s)][k] for s in NC + D1 for k in m6actrl if k in data[("m6a", s)]]
nd1, tt1 = pool(v_rec_d1); nd0, tt0 = pool(v_rec_nc)
ndb, ttb = pool(bg_rec);      ndm, ttm = pool(m6_rec)
w(f"Fisher verdict-Day1 vs inosine-bg: {fisher_exact([[nd1,tt1-nd1],[ndb,ttb-ndb]])[1]:.2e}")
w(f"Fisher verdict-NC   vs inosine-bg: {fisher_exact([[nd0,tt0-nd0],[ndb,ttb-ndb]])[1]:.2e}")
w(f"Fisher verdict-Day1 vs m6A-ctrl:  {fisher_exact([[nd1,tt1-nd1],[ndm,ttm-ndm]])[1]:.2e}")
w(f"Fisher verdict-NC   vs verdict-Day1: {fisher_exact([[nd0,tt0-nd0],[nd1,tt1-nd1]])[1]:.2e}")
sv = [share(r) for r in v_rec_d1 if share(r) is not None]
sb = [share(r) for r in bg_rec if share(r) is not None]
sm = [share(r) for r in m6_rec if share(r) is not None]
w(f"MWU verdict-Day1 vs bg: {mannwhitneyu(sv, sb, alternative='greater')[1]:.2e}")
w(f"MWU verdict-Day1 vs m6A: {mannwhitneyu(sv, sm, alternative='greater')[1]:.2e}")

with open(f"{MM}/mismatch_report.txt", "w") as out:
    out.write("\n".join(report) + "\n")
