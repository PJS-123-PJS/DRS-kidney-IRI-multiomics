#!/usr/bin/env python3
# =============================================================================
# GSE269622 Visium six-timepoint projection (Sham/4h/12h/2d/14d/42d) of the
# seven validated direction sets — extends current Fig 3F (3 timepoints).
# Each timepoint = 1 section (descriptive trajectory, no replicates).
# PT compartment defined by PT-marker module score (top quartile per section).
# Validation anchor: Sham/12h/14d PT m6A_UP trend should reproduce the
# canonical +0.445 / +1.156 / -0.060 pattern from the server pipeline.
# Output: fig3_8mod_spatial/data/GSE269622_visium_6tp_projection_*.csv
# =============================================================================
import pandas as pd, numpy as np, scipy.io as sio, scipy.sparse as sp
import tarfile, gzip, io, os, re, sys

BASE = "/mnt/f/spatial_ext/GSE269622"
GS   = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/upload_5mod_validated"
OUT  = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"
RAW  = f"{BASE}/GSE269622_RAW.tar"

SAMPLES = {  # inner tar member -> timepoint label
    "GSM8323120_visium_sham_male_spaceranger.tar.gz":  "Sham",
    "GSM8323121_visium_hour4_male_spaceranger.tar.gz": "IRI-4h",
    "GSM8323122_visium_hour12_male_spaceranger.tar.gz":"IRI-12h",
    "GSM8323123_visium_day2_male_spaceranger.tar.gz":  "IRI-2d",
    "GSM8323124_visium_day14_male_spaceranger.tar.gz": "IRI-14d",
    "GSM8323125_visium_week6_male_spaceranger.tar.gz": "IRI-42d",
}
TP_ORDER = ["Sham","IRI-4h","IRI-12h","IRI-2d","IRI-14d","IRI-42d"]

SETS = {"m6A_UP":"M6A_UP.txt","m6A_DOWN_res":"M6A_DOWN.txt","m5C_DOWN":"M5C_DOWN.txt",
        "PSU_UP":"PSU_UP.txt","INOSINE_UP":"INOSINE_UP.txt",
        "NM_UP":"NM_UP.txt","NM_DOWN":"NM_DOWN.txt","UM_UP":"UM_UP.txt"}

PT_MARKERS = ["Slc34a1","Lrp2","Ggt1","Slc13a3","Slc22a6","Slc22a8","Cyp7b1",
              "Kap","Aldob","Ass1","Slc5a2","Slc7a13","Acy3","Cndp2","G6pc"]

def read_stream(tarf, member_name, binary=False):
    f = tarf.extractfile(member_name)
    if f is None: raise KeyError(member_name)
    return f.read() if binary else io.TextIOWrapper(f, encoding="utf-8", errors="replace")

def find_members(sample_tar):
    """locate filtered matrix trio inside one spaceranger tar.gz"""
    names = sample_tar.getnames()
    trio = {}
    for n in names:
        ln = n.lower()
        if "filtered_feature_bc_matrix" in ln:
            if ln.endswith("matrix.mtx.gz"):   trio["mtx"] = n
            elif ln.endswith("barcodes.tsv.gz"): trio["bc"] = n
            elif re.search(r"(features|genes)\.tsv\.gz$", ln): trio["ft"] = n
    return trio

results = []
raw = tarfile.open(RAW, "r:")
for member, tp in SAMPLES.items():
    print(f"--- {tp} ({member})")
    inner_bytes = raw.extractfile(member).read()
    st = tarfile.open(fileobj=io.BytesIO(inner_bytes), mode="r:gz")
    trio = find_members(st)
    if len(trio) < 3: print("  !! trio not found:", trio); continue
    ft  = pd.read_csv(read_stream(st, trio["ft"]), sep="\t", header=None,
                      compression="gzip").iloc[:, -1].str.strip()
    bc  = pd.read_csv(read_stream(st, trio["bc"]), sep="\t", header=None,
                      compression="gzip")[0].str.strip()
    mtx = sio.mmread(io.BytesIO(read_stream(st, trio["mtx"], binary=True)))
    M = sp.csr_matrix(mtx)
    if M.shape[0] == len(ft): M = M.T.tocsr()          # to spots x genes
    # log-normalize (CP10k + log1p) — no SCT dependency
    lib = np.asarray(M.sum(axis=1)).ravel(); lib[lib == 0] = 1
    X = M.multiply(1e4 / lib[:, None]).tocsr(); X.data = np.log1p(X.data)
    gidx = pd.Series(np.arange(len(ft)), index=ft.values)
    gmean = np.asarray(X.mean(axis=1)).ravel()
    def mscore(genes):
        genes = [g for g in genes if g in gidx.index]
        if len(genes) < 3: return np.full(X.shape[0], np.nan), 0
        return np.asarray(X[:, gidx.loc[genes].values].mean(axis=1)).ravel() - gmean, len(genes)
    pt_score, npt = mscore(PT_MARKERS)
    thr = np.nanquantile(pt_score, 0.75)
    is_pt = pt_score >= thr
    row = {"timepoint": tp, "n_spots": X.shape[0], "n_PTspots": int(is_pt.sum())}
    for sname, fn in SETS.items():
        genes = pd.read_csv(f"{GS}/{fn}", header=None, sep="\t")[0].str.strip()
        sc, ng = mscore(list(genes))
        z = (sc - np.nanmean(sc)) / np.nanstd(sc)
        row[f"{sname}_n"] = ng
        row[f"{sname}_PT_Z"] = np.nanmean(z[is_pt])
        row[f"{sname}_all_Z"] = np.nanmean(z)
    results.append(row)
    print(f"   spots={row['n_spots']} PTspots={row['n_PTspots']} m6A_UP PT_Z={row['m6A_UP_PT_Z']:+.3f}")

df = pd.DataFrame(results)
df["timepoint"] = pd.Categorical(df["timepoint"], TP_ORDER, ordered=True)
df = df.sort_values("timepoint")
df.to_csv(f"{OUT}/GSE269622_visium_6tp_projection_PT.csv", index=False)
pz = df[["timepoint"] + [c for c in df.columns if c.endswith("_PT_Z")]]
print("\n=== PT-compartment Z across six timepoints ===")
print(pz.to_string(index=False))
print("\nwritten: GSE269622_visium_6tp_projection_PT.csv")
