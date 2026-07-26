# Result file inventory

This inventory lists files that directly support the manuscript claims and the cleaned public reproducibility package.

## Public-package supplied deliverables

- `deliverables/supplementary_tables_v1.xlsx` — consolidated supplementary table workbook. (Handled as a journal supplementary file; not stored in this repository.)

## Core manuscript evidence

- `results/bulk_validation/bulk_validation_enhanced_stats.csv` — final cross-cohort effect sizes, Mann-Whitney p values and FDR q values.
- `results/bulk_validation/bulk_validation_enhanced_module_summary.csv` — module-level replication summary.
- `results/recovery/geo_10x_full_siteaware_gate_decision.json` — TAURUS whole-biopsy recovery gate decision.
- `results/recovery/geo_10x_full_siteaware_recovery_synchrony_tests.csv` — TAURUS recovery tests.

## Bulk cohort outputs

- `results/bulk_validation/GSE16879_bulk_module_scores_clean.csv`
- `results/bulk_validation/GSE16879_bulk_module_tests.csv`
- `results/bulk_validation/GSE23597_bulk_module_scores_clean.csv`
- `results/bulk_validation/GSE23597_bulk_module_tests.csv`
- `results/bulk_validation/GSE14580_bulk_module_scores_clean.csv`
- `results/bulk_validation/GSE14580_bulk_module_tests.csv`

## TAURUS feasibility outputs

- `results/feasibility/GSE282122_sample_metadata.csv`
- `results/feasibility/GSE282122_umap_sample_coverage.csv`
- `results/feasibility/GSE282122_umap_sample_coverage_summary.json`
- `results/feasibility/GSE282122_processed_tar_audit.json`
- `results/feasibility/GSE282122_processed_tar_members.csv`
- `results/module_scores/geo_10x_full_sample_module_scores.csv`
- `results/module_scores/geo_10x_full_10x_module_score_audit.json`

## TAURUS lineage-resolved outputs (2026-07-26 update)

- `results/taurus_lineage/myeloid_sample_state_module_scores.csv`
- `results/taurus_lineage/myeloid_lineage_module_by_state.csv`
- `results/taurus_lineage/epicolonic_patient_state_module_scores.csv`
- `results/taurus_lineage/epicolonic_lineage_module_by_state.csv`
- `results/taurus_lineage/fibperi_patient_state_module_scores.csv`
- `results/taurus_lineage/fibperi_lineage_module_by_state.csv`
- `results/taurus_lineage/composition_within_lineage_proportions.csv`
- `results/taurus_lineage/composition_infiltration_proxies.csv`
- `results/taurus_lineage/composition_infiltration_proxies_colon_only.csv`

## Figures

- `results/figures/fig1_public_data_workflow_schematic.*` — workflow schematic.
- `results/figures/fig_bulk_GSE16879_pretreatment_modules.*` — GSE16879 module figure.
- `results/figures/fig_bulk_GSE23597_baseline_modules.*` — GSE23597 module figure.
- `results/figures/fig_bulk_validation_forest.*` — cross-cohort effect-size summary.
- `results/figures/fig1_geo_10x_full_siteaware_recovery_synchrony.*` — TAURUS recovery-gate figure.
- `results/figures/figS1_feasibility_coverage.*` — sample/cell coverage audit.
- `results/figures/fig5_taurus_composition_key_hits.png` — TAURUS composition key hits (Figure 5 candidate).
- `results/figures/figS3_taurus_lineage_cliffs_delta_heatmap.png` — lineage Cliff's delta heatmap (Supplementary Figure 3 candidate).

Binary figure exports are kept out of this lightweight repository; regenerate them with the `scripts/plot_*.py` scripts from the derived tables.

## Documentation

- `docs/reproducibility_readme.md`
- `docs/source_references.md`
- `docs/analysis_workflow.md`
- `docs/result_file_inventory.md`
- `docs/taurus_zenodo_acquisition_2026-07-26.md`

## Build and reproducibility scripts

- `scripts/parse_gse282122_soft.py`
- `scripts/audit_umap_sample_coverage.py`
- `scripts/audit_gse282122_processed_tar.py`
- `scripts/score_10x_h5_modules_from_tar.py`
- `scripts/compute_recovery_synchrony.py`
- `scripts/summarize_recovery_gate.py`
- `scripts/score_bulk_series_matrix_modules.py`
- `scripts/analyze_gse16879_bulk_modules.py`
- `scripts/analyze_gse23597_bulk_modules.py`
- `scripts/summarize_bulk_validation.py`
- `scripts/enhance_bulk_validation_stats.py`
- `scripts/audit_bulk_module_score_sensitivity.py`
- `scripts/plot_workflow_schematic.py`
- `scripts/plot_gse16879_bulk_validation.py`
- `scripts/plot_gse23597_bulk_validation.py`
- `scripts/plot_bulk_validation_forest.py`
- `scripts/plot_recovery_synchrony.py`
- `scripts/plot_feasibility_coverage.py`
- `scripts/taurus_lineage_module_localization.py`
- `scripts/taurus_cell_composition_by_outcome.py`

Submission-oriented DOCX generation, final packaging and journal-administration scripts are intentionally excluded from the public repository package.
