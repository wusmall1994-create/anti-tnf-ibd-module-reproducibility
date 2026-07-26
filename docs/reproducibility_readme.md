# Reproducibility README

Project: public-database mucosal module analysis of anti-TNF outcomes in inflammatory bowel disease

Last updated: 2026-07-26

> Update 2026-07-26: annotated TAURUS lineage H5AD objects (myeloid, fibroblast/pericyte, epithelial/colonic) were obtained from Zenodo record 10.5281/zenodo.14007626 (CC BY 4.0, MD5-verified; see `docs/taurus_zenodo_acquisition_2026-07-26.md`). Lineage-resolved module localization and cell-composition analyses are now part of the package (`scripts/taurus_lineage_module_localization.py`, `scripts/taurus_cell_composition_by_outcome.py`, outputs under `results/taurus_lineage/`). The GEO whole-biopsy track below remains valid as the feasibility gate; the earlier limitation to GEO matrices alone no longer applies to the lineage-resolved analyses.

## What this repository/package contains

This project contains scientific-analysis scripts, derived tables, figures and reproducibility documentation for a public-database analysis of pretreatment mucosal modules associated with anti-TNF outcomes in IBD.

The analysis has three major evidence tracks:

1. TAURUS/GSE282122 feasibility gate: verifies public single-cell data availability, supports whole-biopsy module scoring and defines the boundary of analyses based on GEO processed matrices.
2. Bulk validation: tests prespecified mucosal modules across GSE16879 and GSE23597, with GSE14580 retained as an overlap audit rather than independent validation.
3. TAURUS lineage-resolved localization (2026-07-26 update): per-cell module scoring and cell-state composition analysis within annotated myeloid, epithelial/colonic and fibroblast/pericyte lineage objects from Zenodo.

## Recommended environment

The analysis was run locally on Windows using Python 3.12.13.

Observed package versions:

- pandas 3.0.3
- numpy 2.5.1
- scipy 1.18.0
- h5py 3.16.0
- matplotlib 3.11.0
- seaborn 0.13.2
- statsmodels 0.14.6
- scikit-learn 1.9.0
- anndata 0.13.1

Install package requirements from `requirements_reproducibility.txt`. A local virtual environment such as `.venv` can be used, but it is not included in the release.

The scientific analysis and figure scripts were smoke-tested from a temporary copy of `public_repository_ready_v2` on 2026-07-19 using the local project virtual environment. The test reran `scripts/summarize_bulk_validation.py`, `scripts/plot_bulk_validation_forest.py` and `scripts/plot_workflow_schematic.py`. Document-generation utilities may require additional packages such as `python-docx`; they are not required for reproducing the scientific analysis tables and figures.

## External data locations

Large input data are stored outside the project directory and must be downloaded from public repositories. For the TAURUS processed GEO archive, use a local path of your choice, for example:

- `<LOCAL_DATA_DIR>/GSE282122_filtered_processed_data.tar.gz`

Public data sources:

- GSE282122: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282122
- TAURUS Zenodo record: https://zenodo.org/records/14007626
- GSE16879: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16879
- GSE23597: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23597
- GSE14580: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14580

## Reproduction order

### 1. TAURUS metadata and archive feasibility

Parse and audit TAURUS metadata:

```powershell
python scripts\parse_gse282122_soft.py
python scripts\audit_umap_sample_coverage.py
python scripts\audit_gse282122_processed_tar.py
```

Key outputs:

- `results/feasibility/GSE282122_sample_metadata.csv`
- `results/feasibility/GSE282122_umap_sample_coverage.csv`
- `results/feasibility/GSE282122_processed_tar_audit.json`
- `results/feasibility/feasibility_gate_report.md`

### 2. TAURUS whole-biopsy module scoring and recovery gate

Run full 10x scoring and site-aware recovery:

```powershell
powershell -ExecutionPolicy Bypass -File scripts\run_geo_10x_full_gate.ps1
```

Main underlying scripts:

```powershell
python scripts\score_10x_h5_modules_from_tar.py
python scripts\compute_recovery_synchrony.py --pair-by-site
python scripts\summarize_recovery_gate.py
python scripts\plot_recovery_synchrony.py
```

Key outputs:

- `results/module_scores/geo_10x_full_sample_module_scores.csv`
- `results/recovery/geo_10x_full_siteaware_patient_lineage_module_recovery.csv`
- `results/recovery/geo_10x_full_siteaware_recovery_synchrony_tests.csv`
- `results/recovery/geo_10x_full_siteaware_gate_decision.json`
- `results/figures/fig1_geo_10x_full_siteaware_recovery_synchrony.*`

Interpretation: whole-biopsy recovery was computable but weak and not suitable as the primary positive result.

### 3. Bulk GEO validation

Download/process GEO series matrices as needed:

```powershell
python scripts\download_geo_series_matrix.py
```

Score and analyse bulk cohorts:

```powershell
python scripts\score_bulk_series_matrix_modules.py
python scripts\analyze_gse16879_bulk_modules.py
python scripts\analyze_gse23597_bulk_modules.py
python scripts\summarize_bulk_validation.py
python scripts\enhance_bulk_validation_stats.py
```

Key outputs:

- `results/bulk_validation/GSE16879_bulk_module_scores_clean.csv`
- `results/bulk_validation/GSE16879_bulk_module_tests.csv`
- `results/bulk_validation/GSE23597_bulk_module_scores_clean.csv`
- `results/bulk_validation/GSE23597_bulk_module_tests.csv`
- `results/bulk_validation/GSE14580_bulk_module_tests.csv`
- `results/bulk_validation/bulk_validation_enhanced_stats.csv`
- `results/bulk_validation/bulk_validation_enhanced_module_summary.csv`

Interpretation: GSE16879 and GSE23597 support the pretreatment module-response association. GSE14580 overlaps with GSE16879 UC results and should not be counted as an independent validation cohort.

### 4. TAURUS lineage-resolved analyses (2026-07-26 update)

Download the annotated lineage H5AD objects from Zenodo record 10.5281/zenodo.14007626 into a local data directory, then run:

```powershell
python scripts\taurus_lineage_module_localization.py
python scripts\taurus_cell_composition_by_outcome.py
```

Key outputs under `results/taurus_lineage/`:

- `myeloid_sample_state_module_scores.csv`, `myeloid_lineage_module_by_state.csv`
- `epicolonic_patient_state_module_scores.csv`, `epicolonic_lineage_module_by_state.csv`
- `fibperi_patient_state_module_scores.csv`, `fibperi_lineage_module_by_state.csv`
- `composition_within_lineage_proportions.csv`, `composition_infiltration_proxies.csv`

Interpretation: within-lineage module expression did not differ significantly by outcome after FDR control; CD remission patients showed a significantly higher within-myeloid proportion of C1Qhi IL1Blo macrophages (p = 6.7e-4, q = 0.038). Cross-lineage proportion comparisons in CD were not interpretable because CD sampling was predominantly ileal.

### 5. Figure generation

Generate manuscript figures:

```powershell
python scripts\plot_workflow_schematic.py
python scripts\plot_gse16879_bulk_validation.py
python scripts\plot_gse23597_bulk_validation.py
python scripts\plot_bulk_validation_forest.py
python scripts\plot_recovery_synchrony.py
python scripts\plot_feasibility_coverage.py
```

Key outputs:

- `results/figures/fig1_public_data_workflow_schematic.*`
- `results/figures/fig_bulk_GSE16879_pretreatment_modules.*`
- `results/figures/fig_bulk_GSE23597_baseline_modules.*`
- `results/figures/fig_bulk_validation_forest.*`
- `results/figures/fig1_geo_10x_full_siteaware_recovery_synchrony.*`
- `results/figures/figS1_feasibility_coverage.*`
- `results/figures/fig5_taurus_composition_key_hits.png`
- `results/figures/figS3_taurus_lineage_cliffs_delta_heatmap.png`

Binary figure exports are kept out of this lightweight repository; the scripts above regenerate them from the derived tables.

### 6. Supplementary workbook

The supplementary table workbook is handled as a journal supplementary file and is not stored in this repository.

Submission-oriented Word document generation and ZIP finalization scripts are intentionally not included in the public repository package. They are administrative packaging tools rather than scientific-analysis dependencies.

## Primary result files for manuscript claims

- Main cross-cohort statistics: `results/bulk_validation/bulk_validation_enhanced_stats.csv`
- Module summary: `results/bulk_validation/bulk_validation_enhanced_module_summary.csv`
- TAURUS feasibility-gate decision: `results/recovery/geo_10x_full_siteaware_gate_decision.json`
- TAURUS lineage composition statistics: `results/taurus_lineage/composition_within_lineage_proportions.csv`

## Claims supported by current outputs

Supported:

- Pretreatment mucosal module activity is associated with later anti-TNF response in public IBD cohorts.
- Myeloid inflammation is the most consistent module, lower in responders across all independent comparisons and meeting FDR q < 0.10 in all tested strata.
- TAURUS whole-biopsy module scoring is feasible, but whole-biopsy recovery toward healthy is not robust enough as a primary response metric.
- GSE14580 should be treated as an overlap audit rather than an independent validation cohort.
- Within TAURUS myeloid lineage, baseline composition differs by outcome in CD (higher C1Qhi IL1Blo macrophage proportion in remission patients).

Not supported without further data:

- A clinically deployable diagnostic test.
- A new causal anti-TNF resistance mechanism.
- Within-lineage module expression differences by outcome in TAURUS (no comparison survived FDR control).
- Cross-lineage single-cell synchrony at the patient level.

Do not upload very large raw archives unless the repository is designed for large data. Instead, cite GEO/Zenodo accessions and document the expected external data path with a placeholder such as `<LOCAL_DATA_DIR>`.
