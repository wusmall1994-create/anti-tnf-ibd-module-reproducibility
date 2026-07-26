# TAURUS Zenodo annotated H5AD acquisition record

Date: 2026-07-26
Source: https://doi.org/10.5281/zenodo.14007626 (record 14007626, latest version, 2024-10-30)
License: CC BY 4.0 (open access)
Download method: curl with --resolve DNS bypass (local DNS failed to resolve zenodo.org; public DNS OK)
Destination: data/taurus_zenodo/

| File | Bytes | MD5 (Zenodo record) | MD5 verified |
|---|---:|---|---|
| myeloid_final.h5ad | 416,722,961 | bdfe50345a11abdb1a72b2439bf9950e | PASS |
| fibperi_final.h5ad | 1,407,900,848 | 62bfe745f6c7892d7c1069e220136bfc | PASS |
| epicolonic_final.h5ad | 6,031,357,405 | 99072da25c47bc797ceba8d8bd5d8c5a | PASS |
| paired_sample_list.csv | 3,102 | 3300a53889bb4b70c48ec66dbb66beea | size match |

Annotation structure (verified via anndata backed read, anndata 0.13.2):
- myeloid: 30,858 cells x 33,075 genes; 11 final_analysis cell states (e.g. S100A8 A9hi mono, C1Qhi IL1Bhi macro, LAMP3pos DC)
- fibperi: 103,053 cells x 33,075 genes
- epicolonic: 305,862 cells x 33,075 genes
- obs columns include: sample_id, Patient, Disease, Site, Treatment (Pre/Post), Inflammation_score, Remission_status, final_analysis/minor/major/sub_bucket/bucket
- Zenodo record note: for baseline analyses, use baseline samples with inflammation score > 6.5 to recapitulate paper analyses

Purpose: lineage-level localization of the four prespecified modules (myeloid inflammation, epithelial IFN/damage, epithelial barrier maturity, stromal fibroinflammatory) to upgrade the manuscript's TAURUS section from a feasibility gate to true cell-state analysis.
