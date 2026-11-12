# DRS-kidney-IRI-multiomics

Analysis code accompanying the manuscript:

> **Single-molecule stoichiometry by direct RNA sequencing distinguishes epitranscriptome remodeling from artifacts in kidney ischemia–reperfusion injury**
> (under review)

This repository contains the figure-generation and analysis scripts (R / Python / Bash) used to produce every main-figure and supplementary-figure result of the study — a five-group time course (Sham / IS / Day 1 / Day 3 / Day 7, n = 3 mice per group) of mouse kidney ischemia–reperfusion injury (IRI) profiled by Oxford Nanopore direct RNA sequencing (DRS), matched short-read RNA-seq, ATAC-seq and untargeted metabolomics.

## Data availability (companion datasets)

| Data type | Repository | Accession |
|:---|:---|:---|
| DRS + RNA-seq + ATAC-seq reads (75 files) | NGDC **GSA** (BioProject [PRJCA079331](https://ngdc.cncb.ac.cn/bioproject/browse/PRJCA079331), submission subCRA087106) | GSA accession pending assignment |
| Untargeted metabolomics (raw mzXML + results) | NGDC **OMIX** | [OMIX023403](https://ngdc.cncb.ac.cn/omix/) |
| Single-cell RNA-seq (this study) | GEO | [GSE274819](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274819) |
| Public re-analysis cohorts | GEO | GSE139107, GSE171417, GSE252106, GSE151302, GSE269622, GSE298899, GSE224289 |

## Repository layout

Each directory holds the scripts of one figure; `MANIFEST.csv` at the root maps every script to its figure/panel(s) with MD5 checksums.

| Directory | Content |
|:---|:---|
| `fig1_platform_qc/scripts/` | Fig 1 — platform overview & data quality control |
| `fig2_8mod_temporal/scripts/` | Fig 2 — temporal dynamics of the multi-modification landscape |
| `fig3_8mod_spatial/scripts/` | Fig 3 — stoichiometry precision & spatial/cell-type projection, public-data orthogonal validation |
| `fig4_as_landscape/scripts/` | Fig 4 — long-read alternative-splicing landscape |
| `fig5_polya_pt/scripts/` | Fig 5 — poly(A)-tail & proximal-tubule cell-type integration |
| `fig6_novel_transcript/scripts/` | Fig 6 — novel transcripts / dark transcriptome |
| `fig7_multiomics/scripts/` | Fig 7 — multi-omics integration & temporal decoupling |
| `fig8_validation/scripts/` | Fig 8 — experimental validation |
| `supplement/scripts/` | Supplementary figures/tables build scripts |
| `shared/` | Shared palettes and utilities |

## Requirements & usage

- **R** (≥ 4.2) with the packages loaded at the top of each script (tidyverse/ggplot2 stack plus Bioconductor packages where indicated); **Python 3** for the mismatch-extraction and spatial-projection utilities; Bash drivers assume a Linux environment.
- Scripts were run against local copies of the deposited data and reference absolute working paths (e.g. `E:/`, `F:/`). After downloading the data from GSA/OMIX/GEO, adapt the path variables at the top of each script to your local layout.
- Raw sequencing reads are intentionally **not** part of this repository — they live in GSA/OMIX under the accessions above.

## Authors

Jiashan Pan, Mingchen Cai, Qi Zhu, Mo Yang, Xinlei Ling, Wei Chen, Handong Ding, Jinbiao Zhong, Zongyao Hao\*, Chaozhao Liang\*, Xiaoming Meng\*, Guiyi Liao\* (\*corresponding authors)

The First Affiliated Hospital of Anhui Medical University & Anhui Medical University, Hefei, China

## Contact & citation

Questions: Jiashan Pan (panjs@fy.ahmu.edu.cn) or the corresponding authors (liaoguiyi@ahmu.edu.cn).

If you use this code, please cite the manuscript and the Zenodo archive of this repository (DOI minted on first release).

## License

MIT — see [LICENSE](LICENSE).
