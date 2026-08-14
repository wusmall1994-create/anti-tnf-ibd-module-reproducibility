# Cross-cohort calibration of pretreatment mucosal modules in IBD

This repository contains the analysis code, reproducibility documentation and
selected derived outputs supporting a public-data transcriptomic study of
pretreatment mucosal states associated with biologic response in inflammatory
bowel disease.

Version 2.0.0 adds dataset-independence auditing, a unified benchmark of eight
published signatures, leave-one-cohort-out and stratified analyses, clinical
increment modelling, TAURUS-referenced partial deconvolution and treatment-
boundary analyses beyond anti-TNF therapy.

## Repository contents

- `scripts/`: data acquisition, harmonisation, scoring, meta-analysis,
  deconvolution and figure-generation code.
- `results/`: selected machine-readable derived CSV and JSON tables, source
  data and figure exports.
- `docs/`: workflow, provenance, reproducibility and release documentation.
- `requirements_reproducibility.txt`: Python environment requirements.
- `CITATION.cff`, `LICENSE` and `zenodo_metadata_template.json`: citation,
  licensing and deposition metadata.

The journal manuscript, cover letter, submission forms, correspondence and
local packaging files are intentionally excluded.

## Public datasets

The analyses use publicly available data from:

- GSE282122 and Zenodo 10.5281/zenodo.14007626: TAURUS single-cell atlas and
  annotated lineage objects.
- GSE16879, GSE23597, GSE92415 and E-MTAB-7604: primary anti-TNF bulk cohorts.
- GSE14580 and GSE12251: overlap and independence auditing.
- GSE73661: vedolizumab treatment-boundary analysis.
- GSE206285: ustekinumab treatment-boundary analysis.

Raw public archives are not redistributed. Accession links and acquisition
details are listed in `docs/source_references.md` and the workflow documents.

## Main v2.0.0 entry points

- `scripts/run_calibration_upgrade.py`
- `scripts/run_taurus_reference_deconvolution.py`
- `results/calibration_upgrade_20260814/`
- `docs/release_notes_v2.0.0.md`

The upgrade pipeline covers:

1. accession- and expression-level independence auditing;
2. reconstruction and external benchmarking of eight published signatures;
3. leave-one-cohort-out and disease, tissue and drug strata;
4. repeated cross-validated GSE92415 clinical-increment models;
5. TAURUS-referenced five-state partial deconvolution and sensitivity analysis;
6. vedolizumab and ustekinumab treatment-boundary analysis.

## Environment and local data layout

The workflow was run with Python 3.12.13. Install the documented dependencies:

```bash
pip install -r requirements_reproducibility.txt
```

Place downloaded public datasets under a local analysis root following the
paths documented in `docs/reproducibility_readme.md`. Set the root before
running the v2.0.0 upgrade scripts:

```bash
export ANTI_TNF_IBD_ROOT=/path/to/analysis-root
python scripts/run_calibration_upgrade.py
python scripts/run_taurus_reference_deconvolution.py
```

On PowerShell:

```powershell
$env:ANTI_TNF_IBD_ROOT = '<analysis-root>'
python scripts/run_calibration_upgrade.py
python scripts/run_taurus_reference_deconvolution.py
```

If `ANTI_TNF_IBD_ROOT` is not set, the repository root is used.

## Privacy and provenance

This release contains no directly identifiable participant data, manuscript
files, submission correspondence, email addresses, telephone numbers, local
absolute paths, credentials or private keys. Public de-identified sample
identifiers are retained where needed for reproducibility. Public author names
and institutional affiliations appear only in citation and repository metadata.

## License and citation

Code is released under the MIT License. Cite the software using `CITATION.cff`
and cite the original public datasets and publications listed in the repository
documentation.
