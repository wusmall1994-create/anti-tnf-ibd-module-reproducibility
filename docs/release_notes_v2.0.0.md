# Release notes: v2.0.0

This release expands the original reproducibility archive into a cohort-audited
cross-study calibration resource. It contains code, documentation, selected
derived tables and figure outputs. It does not contain manuscript files,
submission correspondence, private clinical data or local administrative files.

## New analyses

1. Dataset-independence audit using accession metadata, sample identifiers,
   endpoint definitions, SHA-256 file hashes and expression fingerprints.
2. Unified reconstruction and external benchmarking of eight published
   mucosal signatures under a common preprocessing and meta-analysis pipeline.
3. Leave-one-cohort-out analyses and prespecified disease, tissue and drug
   strata.
4. Repeated cross-validated clinical-increment analysis in GSE92415, including
   AUC, Brier score, calibration, NRI, IDI and decision-curve outputs.
5. TAURUS-referenced five-state partial-reference deconvolution, marker-set
   sensitivity analysis, synthetic-mixture recovery and composition-adjusted
   modelling.
6. Treatment-boundary analyses in vedolizumab-treated GSE73661 and
   ustekinumab-treated GSE206285.

## Main entry points

- `scripts/run_calibration_upgrade.py`
- `scripts/run_taurus_reference_deconvolution.py`
- `results/calibration_upgrade_20260814/`
- `docs/analysis_workflow.md`
- `docs/reproducibility_readme.md`

## Local data paths

Raw public datasets are not redistributed. Set `ANTI_TNF_IBD_ROOT` to the local
analysis root containing the documented `data/` layout before running the two
upgrade scripts. If the variable is not set, the repository root is used.

## Privacy and scope

The release was screened before publication for manuscript and office-document
files, correspondence, email addresses, telephone numbers, absolute local paths,
credentials and private keys. Public author names and institutional affiliations
are retained only in citation and repository metadata.
