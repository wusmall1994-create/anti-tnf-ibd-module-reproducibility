# Interpretable anti-TNF IBD mucosal module analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21676792.svg)](https://doi.org/10.5281/zenodo.21676792)

This repository contains cleaned code, selected derived result tables, and figure exports for a public-database transcriptomic analysis of pretreatment mucosal modules associated with anti-TNF outcomes in inflammatory bowel disease.

It is a scientific reproducibility repository, not the journal-administration package. Manuscript Word files, cover letters, submission forms, response letters, and local packaging/audit administration files are intentionally excluded.

## Contents

This repository includes:

- analysis scripts under `scripts/`
- reproducibility documentation under `docs/`
- selected derived CSV/JSON result tables under `results/`
- figure-generation scripts under `scripts/` (`plot_*.py` regenerate all manuscript PNG/SVG figures from the derived tables; binary figure exports are kept out of this lightweight repository)
- Python package requirements in `requirements_reproducibility.txt`
- citation/license/deposition metadata templates

The supplementary Excel workbook is handled as a journal supplementary file rather than stored in this GitHub repository, because the repository is kept focused on code, machine-readable derived tables, and figures.

## Public datasets

- GSE282122: TAURUS anti-TNF single-cell atlas.
- GSE16879: mucosal expression profiling before and after first infliximab treatment.
- GSE23597: colonic biopsy expression data from infliximab-treated ulcerative colitis patients.
- GSE14580: ulcerative colitis infliximab response dataset retained as an overlap audit.

The large TAURUS GEO processed archive is not included. It should be downloaded from public GEO accession GSE282122:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282122

When rerunning the TAURUS 10x scoring step, place the downloaded archive in a local data directory of your choice and pass that path to the relevant script or PowerShell wrapper.

## Quick reproduction guide

See:

- `docs/reproducibility_readme.md`
- `docs/result_file_inventory.md`

Core manuscript-supporting outputs include:

- `results/bulk_validation/bulk_validation_enhanced_stats.csv`
- `results/bulk_validation/bulk_validation_enhanced_module_summary.csv`
- `results/bulk_validation/bulk_module_score_sensitivity_summary.json`
- `results/recovery/geo_10x_full_siteaware_gate_decision.json`

## Environment

The workflow was run using Python 3.12.13. Install packages from:

```bash
pip install -r requirements_reproducibility.txt
```

## Citation

The exact software release used for the manuscript is archived in Zenodo:

> Liang H, Wu X, Su X, Qiu R, Zheng F. *Interpretable public-data validation of pretreatment mucosal modules associated with anti-TNF outcomes in inflammatory bowel disease*. Version 1.0.0 [software]. Zenodo. 2026. https://doi.org/10.5281/zenodo.21676792

The version-specific DOI should be used when citing the analysis code. Release `v1.0.0` corresponds to commit `49cbf3fdc8436d5f96d3f7c3d6d567f88a26ce96`.

## Scientific boundary

The current release supports a public transcriptomic module-validation study. It does not support claims of:

- a clinically deployable diagnostic test;
- a causal anti-TNF resistance mechanism;
- cell-type-specific localization from the current GEO matrices alone;
- cross-lineage single-cell synchrony without annotated TAURUS H5AD files.

## License

Code is released under the MIT License (see `LICENSE`). Documentation and non-code research materials are intended for CC BY 4.0-style reuse unless otherwise specified by the author team or journal policy.

## Lineage-resolved update (2026-07-26)

Annotated lineage-level TAURUS H5AD objects were obtained from Zenodo record
10.5281/zenodo.14007626 (CC BY 4.0; MD5-verified, see
`docs/taurus_zenodo_acquisition_2026-07-26.md`). This update adds:

- `scripts/taurus_lineage_module_localization.py`: per-cell module scoring within
  the myeloid lineage object, aggregated sample -> patient (patient = statistical unit)
- `scripts/taurus_cell_composition_by_outcome.py`: within-lineage cell-state
  proportions and cross-lineage infiltration proxies by outcome
- derived tables under `results/taurus_lineage/`
- figure scripts regenerate `results/figures/fig5_taurus_composition_key_hits.png` (Figure 5 candidate)
  and `results/figures/figS3_taurus_lineage_cliffs_delta_heatmap.png` (Supplementary Figure 3 candidate)

Baseline lineage analyses use pretreatment samples with inflammation score > 6.5,
following the TAURUS authors' guidance.
