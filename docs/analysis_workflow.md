# Reproducible workflow

This document tracks the actual workflow for the public-database IBD anti-TNF project.This document tracks the actual workflow for the public-database IBD anti-TNF project.

## 0. Environment

Use a local project virtual environment:Use a local project virtual environment:

```powershellpowershellpowershellpowershell
.\.venv\Scripts\python.exe
```

Primary packages:Primary packages:

- pandas pandas pandas pandas
- numpy numpy numpy numpy
- scipy scipy scipy scipy
- h5py h5py h5py h5py
- anndata anndata anndata anndata

## 1. Metadata gate

```powershellpowershellpowershellpowershell
.\.venv\Scripts\python.exe scripts\parse_gse282122_soft.py
```

Outputs:Outputs:

- `results/feasibility/GSE282122_sample_metadata.csv`

## 2. Outcome gate

Outcome labels are obtained from:Outcome labels are obtained from:

- `data/metadata/paired_sample_list.csv` (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626) (TAURUS Zenodo record 10.5281/zenodo.14007626)

Current outcome-linked patient counts:Current outcome-linked patient counts:

- UC remission: 4 UC remission: 4 UC remission: 4 UC remission: 4
- UC non-remission: 13 UC non-remission: 13 UC non-remission: 13 UC non-remission: 13
- CD remission: 10 CD remission: 10 CD remission: 10 CD remission: 10
- CD non-remission: 5 CD non-remission: 5 CD non-remission: 5 CD non-remission: 5

## 3. UMAP coverage gate

```powershell
.\.venv\Scripts\python.exe scripts\audit_umap_sample_coverage.py
```

Outputs:

- `results/feasibility/GSE282122_umap_sample_coverage.csv`
- `results/feasibility/GSE282122_umap_sample_coverage_summary.json`

Current result:

- 987,743 cells
- 216 samples
- 0 unmatched samples

## 4. Processed expression data download

Preferred direct NCBI FTP URL:

```text
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282122/suppl/GSE282122_filtered_processed_data.tar.gz
```

Local target:

```text
<LOCAL_DATA_DIR>/GSE282122_filtered_processed_data.tar.gz
```

Check progress:

```powershell
.\.venv\Scripts\python.exe scripts\check_download_status.py
```

## 5. Processed tar audit

Run after the tar file finishes downloading:

```powershell
.\.venv\Scripts\python.exe scripts\audit_gse282122_processed_tar.py
```

Outputs:

- `results/feasibility/GSE282122_processed_tar_audit.json`
- `results/feasibility/GSE282122_processed_tar_members.csv`

## 6. Whole-biopsy 10x module scoring from GEO

The GEO processed archive contains per-sample CellRanger `filtered_feature_bc_matrix.h5` files. These support whole-biopsy / sample-level module scoring but do not contain TAURUS cell-state annotations.

After the tar file finishes downloading:

```powershell
.\.venv\Scripts\python.exe scripts\score_10x_h5_modules_from_tar.py --label geo_10x_full
```

Outputs:

- `results/module_scores/geo_10x_full_sample_module_scores.csv`
- `results/module_scores/geo_10x_full_10x_module_score_audit.json`

This route is the first expression-level gate. It can validate whether public expression data and outcome metadata support a recovery-score analysis at whole-biopsy level.

## 7. Annotated lineage-level module scoring

For each H5AD/AnnData object identified in the processed data:

```powershell
.\.venv\Scripts\python.exe scripts\score_h5ad_modules_by_sample.py --h5ad PATH_TO_OBJECT.h5ad --label myeloid
```

Outputs:

- `results/module_scores/{label}_sample_module_scores.csv`
- `results/module_scores/{label}_sample_celltype_module_scores.csv`
- `results/module_scores/{label}_module_score_audit.json`

This route requires annotated H5AD objects, such as the TAURUS Zenodo lineage objects. It is needed for the full myeloid/epithelial/stromal synchrony claim.

## 8. Recovery and synchrony statistics

After module score CSVs exist:

```powershell
.\.venv\Scripts\python.exe scripts\compute_recovery_synchrony.py `
  --score-csv results\module_scores\myeloid_sample_module_scores.csv `
  --score-csv results\module_scores\epicolonic_sample_module_scores.csv `
  --score-csv results\module_scores\fibperi_sample_module_scores.csv `
  --label first_pass
```

Outputs:

- `results/recovery/first_pass_patient_lineage_module_recovery.csv`
- `results/recovery/first_pass_patient_synchrony.csv`
- `results/recovery/first_pass_recovery_synchrony_tests.csv`

For GEO whole-biopsy scoring, run:

```powershell
.\.venv\Scripts\python.exe scripts\compute_recovery_synchrony.py `
  --score-csv results\module_scores\geo_10x_full_sample_module_scores.csv `
  --label geo_10x_full
```

Note: a single whole-biopsy score file can produce recovery statistics but not cross-lineage synchrony. Synchrony requires at least two or three annotated lineage-level score files.

## 9. TAURUS lineage-resolved module localization (2026-07-26 update)

With annotated lineage H5AD objects from Zenodo record 10.5281/zenodo.14007626:

```powershell
.\.venv\Scripts\python.exe scripts\taurus_lineage_module_localization.py
.\.venv\Scripts\python.exe scripts\taurus_cell_composition_by_outcome.py
```

Outputs under `results/taurus_lineage/`:

- per-lineage sample-level module scores and state-level summaries (`*_lineage_module_by_state.csv`)
- within-lineage cell-state proportions and cross-lineage infiltration proxies by outcome

Baseline analyses use pretreatment samples with inflammation score > 6.5, following the TAURUS authors' guidance. Patients are the statistical unit (scores aggregated sample -> patient).

