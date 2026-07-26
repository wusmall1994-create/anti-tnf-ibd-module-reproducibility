# Public repository package smoke test

Date: 2026-07-19

## Purpose

Verify that the cleaned public repository package can be copied to a fresh temporary directory and used to rerun lightweight scientific-analysis outputs without relying on the manuscript submission package.

## Package tested

- Source directory: `deliverables/public_repository_ready_v2`
- Temporary-copy execution: yes
- Large external TAURUS raw archive required: no for this smoke test

## Environment

The test used the local project virtual environment. In a fresh checkout or archive extraction, create an equivalent environment and install `requirements_reproducibility.txt`. Generic example:

```text
<PROJECT_ROOT>/.venv/Scripts/python.exe
```

Observed scientific-analysis package versions:

```text
pandas==3.0.3
numpy==2.5.1
scipy==1.18.0
h5py==3.16.0
matplotlib==3.11.0
seaborn==0.13.2
statsmodels==0.14.6
scikit-learn==1.9.0
anndata==0.13.1
```

`python-docx` was not present in this virtual environment; it is only needed for generating submission DOCX files and is not required for the scientific-analysis smoke test.

## Commands tested

From a temporary copy of `public_repository_ready_v2`:

```powershell
python scripts\summarize_bulk_validation.py
python scripts\plot_bulk_validation_forest.py
python scripts\plot_workflow_schematic.py
```

## Result

Status: pass.

The summary script printed the expected cohort/module table. The plotting scripts regenerated the following output families in the temporary copy:

- `results/figures/fig_bulk_validation_forest.*`
- `results/figures/fig1_public_data_workflow_schematic.*`

## Interpretation

This smoke test supports that the cleaned public package contains the key derived tables, scripts and figure-generation resources needed for lightweight reproduction of central bulk-validation and workflow outputs. Full TAURUS 10x reruns still require downloading the large public GSE282122 processed archive separately from GEO/Zenodo.
