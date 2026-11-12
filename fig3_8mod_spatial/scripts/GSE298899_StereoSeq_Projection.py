#!/usr/bin/env python3
# =============================================================================
# GSE298899 Stereo-seq (sham n=3 vs bIRI n=3, 14d) projection of the seven
# validated direction sets — orthogonal check for Fig 3F (Visium) extension.
# Author-provided compartment annotation (figure_idents) incl. PT-S1/S2, PT-S3,
# FR-PT. Module score = mean(SCT of set) - mean(SCT of all genes), per bin.
# Z = scale(score) across ALL bins pooled; group means per sample x idents.
# Stats at sample level (n=3/group): mean +- SD, Welch t (direction report).
# Output: fig3_8mod_spatial/data/GSE298899_stereoseq_projection.{tsv,csv}
# =============================================================================
import pandas as pd, numpy as np, scipy.io as sio, scipy.sparse as sp
from scipy.stats import ttest_ind
import os

BASE = "/mnt/f/spatial_ext/GSE298899"
GS   = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data/upload_5mod_validated"
OUT  = "/mnt/d/desktop/DRS/DRS_TechShowcase/fig3_8mod_spatial/data"

SETS = {  # seven validated direction sets (Fig 3F caliber) + key contrasts
    "m6A_UP":     "M6A_UP.txt",
    "m6A_DOWN_res": "M6A_DOWN.txt",
    "m5C_DOWN":   "M5C_DOWN.txt",
    "PSU_UP":     "PSU_UP.txt",
    "INOSINE_UP": "INOSINE_UP.txt",
    "NM_UP":      "NM_UP.txt",
    "NM_DOWN":    "NM_DOWN.txt",
    "UM_UP":      "UM_UP.txt",
}

md  = pd.read_csv(f"{BASE}/GSE298899_metadata.csv.gz")
bc  = pd.read_csv(f"{BASE}/GSE298899_barcodes.tsv.gz", header=None)[0].str.strip()
ft  = pd.read_csv(f"{BASE}/GSE298899_raw_features.tsv.gz", header=None)[0].str.strip()
print(f"bins meta={len(md)} barcodes={len(bc)} features={len(ft)}")

print("reading mtx (this takes a few minutes) ...")
M = sio.mmread(f"{BASE}/GSE298899_sct_matrix.mtx.gz").tocsr()
print("mtx shape:", M.shape)
# orient to bins x genes
if M.shape[0] == len(ft) and M.shape[1] == len(bc):
    M = M.T.tocsr()
elif not (M.shape[0] == len(bc)):
    raise SystemExit("orientation mismatch")

md = md.set_index("barcodes").loc[bc.values].reset_index()
assert (md["barcodes"].values == bc.values).all()

gidx = pd.Series(np.arange(len(ft)), index=ft.values)
global_mean = np.asarray(M.mean(axis=1)).ravel()          # per-bin mean over all genes

scores = {}
for name, fn in SETS.items():
    genes = pd.read_csv(f"{GS}/{fn}", header=None, sep="\t")[0].str.strip()
    genes = [g for g in genes if g in gidx.index]
    sub = M[:, gidx.loc[genes].values]
    scores[name] = np.asarray(sub.mean(axis=1)).ravel() - global_mean
    print(f"  {name}: {len(genes)} genes matched")

S = pd.DataFrame(scores)
Sz = (S - S.mean()) / S.std(ddof=0)                        # pooled Z across all bins
res = pd.concat([md[["barcodes","sample","figure_idents"]], Sz], axis=1)
res["group"] = np.where(res["sample"].str.startswith("sham"), "Sham", "IRI_14d")

# sample x idents x set means
long = res.melt(id_vars=["barcodes","sample","group","figure_idents"],
                var_name="set", value_name="z")
summ = (long.groupby(["set","figure_idents","group","sample"])["z"]
            .mean().reset_index())
pv = []
for (s, comp), gdf in summ.groupby(["set","figure_idents"]):
    a = gdf[gdf.group=="Sham"]["z"].values; b = gdf[gdf.group=="IRI_14d"]["z"].values
    if len(a) >= 2 and len(b) >= 2:
        t, p = ttest_ind(b, a, equal_var=False)
        pv.append((s, comp, a.mean(), b.mean(), b.mean()-a.mean(), p))
pvt = pd.DataFrame(pv, columns=["set","compartment","Sham_meanZ","IRI14d_meanZ",
                                "deltaZ","welch_p"]).sort_values(["set","deltaZ"])
long.to_csv(f"{OUT}/GSE298899_stereoseq_projection_bins.tsv", sep="\t", index=False)
summ.to_csv(f"{OUT}/GSE298899_stereoseq_projection_sample_means.csv", index=False)
pvt.to_csv(f"{OUT}/GSE298899_stereoseq_projection_stats.csv", index=False)

pd.set_option("display.width", 200)
key = ["PT-S1/S2","PT-S3","FR-PT","TAL","CNT","Injured tubule"]
print("\n=== key compartments (Sham vs IRI-14d mean Z) ===")
print(pvt[pvt.compartment.isin(key)].to_string(index=False))
print("\nwritten: GSE298899_stereoseq_projection_{bins.tsv,sample_means.csv,stats.csv}")
