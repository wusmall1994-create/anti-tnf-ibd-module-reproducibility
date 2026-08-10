# Small but substantive enhancement gate

Date: 2026-08-05

## Decision

**PASS, with a controlled claim set.** The enhancement is scientifically substantive enough to justify a manuscript revision. Two additional independent pretreatment mucosal cohorts can be added, formal random-effects meta-analysis is feasible, baseline disease activity can be adjusted in one new cohort, anatomical site can be controlled in the other, and the current modules can be compared quantitatively with published signatures.

This gate does **not** support claiming a new superior predictor or a new single-cell discovery. It supports a stronger cross-cohort validation and compartment-informed synthesis.

## 1. New independent outcome-labelled cohorts

### GSE92415 — include

- PURSUIT-SC golimumab trial, ulcerative colitis, pretreatment colonic mucosa.
- Eligible analysis set: 59 patients (32 responders; 27 non-responders).
- Baseline age and total Mayo score are available.
- The myeloid module had an AUC of 0.721 for higher score predicting non-response and Hedges g of -0.743 (responder minus non-responder).
- Primary record: [NCBI GEO GSE92415](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92415).

### E-MTAB-7604 — include

- Pretreatment inflamed intestinal mucosal RNA-seq, 44 patients (19 responders; 25 non-responders).
- Tissue site is explicit: 27 colon and 17 ileum; drug is explicit: 27 adalimumab and 17 infliximab.
- The myeloid direction reproduced in both colon (g=-0.687) and ileum (g=-0.685).
- Patient-level clinical activity is not publicly available, so tissue and drug were adjusted instead.
- Primary record: [EMBL-EBI BioStudies E-MTAB-7604](https://www.ebi.ac.uk/biostudies/studies/E-MTAB-7604); sequencing project [ENA PRJEB30904](https://www.ebi.ac.uk/ena/browser/view/PRJEB30904).

### Single-cell search

No second truly independent public pretreatment single-cell cohort with anti-TNF outcome labels and a comparable mucosal endpoint was identified. TAURUS should remain the lineage-localisation component, not be portrayed as one of several independent single-cell validations.

### Important exclusions

- GSE12251 is an exact sample-level subset of GSE23597 and cannot be counted independently.
- GSE73661 is conservatively excluded from the independent meta-analysis because its Leuven provenance and 8/15 response composition closely mirror the older Arijs UC cohort; public metadata do not establish patient-level independence.
- Available blood cohorts are biologically independent but are not valid direct tests of mucosal compartment-informed modules.

Full decisions are in `candidate_dataset_audit.csv`; the exact GSE12251/GSE23597 audit is in the existing overlap table.

## 2. Formal meta-analysis

Effect size was Hedges g (responder minus non-responder). Module scores were means of within-stratum gene z-scores. The primary analysis first combined disease/tissue strata within each accession, then fit a four-cohort REML random-effects model with Hartung-Knapp intervals. This prevents the two GSE16879 disease strata and the two E-MTAB tissue strata from being counted as independent studies.

| Module | Pooled g | 95% HK CI | p | I2 |
|---|---:|---:|---:|---:|
| Myeloid inflammation | -0.993 | -1.669 to -0.318 | 0.018 | 33.8% |
| Epithelial barrier maturity | 0.562 | 0.183 to 0.941 | 0.018 | 0.0% |
| Epithelial IFN/damage | -0.623 | -1.536 to 0.289 | 0.118 | 59.7% |
| Stromal fibroinflammatory | -0.630 | -1.409 to 0.150 | 0.082 | 53.7% |

The colon-only sensitivity analysis gave myeloid g=-1.020 (95% HK CI -1.702 to -0.338; p=0.018; I2=28.5%), showing that the result was not created by inclusion of ileal E-MTAB samples.

Interpretation: the myeloid module is the clearest primary cross-cohort result. Barrier maturity is a credible secondary result. The IFN/damage and stromal modules remain directionally suggestive but should not be promoted as statistically established pooled effects.

## 3. Covariate and anatomical control

### GSE92415

After adjustment for baseline total Mayo score and age, the myeloid module remained associated with subsequent response:

- OR for response per 1-SD higher score: 0.474
- 95% CI: 0.244 to 0.921
- p=0.0275

This directly addresses the concern that the module merely recapitulates baseline clinical activity, although residual confounding by histologic inflammation remains possible.

### E-MTAB-7604

After adjustment for ileum/colon site and infliximab/adalimumab treatment:

- Myeloid module: OR for response 0.395 (95% CI 0.173 to 0.903; p=0.0276).
- Barrier maturity module: OR for response 2.206 (95% CI 1.070 to 4.551; p=0.0322).

The tissue-stratified effects and colon-only meta-analysis are preferable to pretending that anatomical heterogeneity is absent.

## 4. Quantitative comparison with published signatures

Benchmarks were prespecified from published work:

- Arijs 2009 top-five genes: TNFRSF11B, STC1, PTGS2, IL13RA2 and IL11.
- West 2017 OSM-associated 22-gene cytokine/chemokine set.
- West 2017 OSM/OSMR two-gene axis.

Signature sources: [Arijs et al., Gut 2009](https://pubmed.ncbi.nlm.nih.gov/19700435/) and [West et al., Nature Medicine 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5420447/).

In GSE92415, the myeloid module AUC was 0.721 versus 0.659, 0.688 and 0.668 for the three benchmarks. Paired stratified bootstrap differences all crossed zero.

In E-MTAB-7604, the myeloid module AUC was 0.714 versus 0.703, 0.648 and 0.739. Again, no benchmark difference was clearly significant; versus West OSM22 the numerical difference was 0.065 (95% bootstrap CI -0.002 to 0.141; p=0.0668).

Therefore the defensible conclusion is **comparable performance with a compact compartment-informed module**, not superiority.

The benchmark comparisons are descriptive rather than clean external validations of every published signature: GSE92415 contributed to the clinical-trial evidence used in the West study, and the older Leuven cohorts contributed to the Arijs evidence. E-MTAB-7604 supplies the most useful additional external comparison, but no public cohort can remove all historical-data reuse concerns.

Gene-set overlap also constrains novelty. The myeloid module shares five genes with the West OSM22 set (CCL3, CCL4, CXCL8, IL1B and OSM; overlap coefficient 0.556). The manuscript should explicitly frame the contribution as integration, robustness and cellular contextualisation rather than a wholly new inflammatory signature.

## Recommended manuscript enhancement

1. Add GSE92415 and E-MTAB-7604 to Methods, dataset table and Results.
2. Make the four-cohort random-effects meta-analysis the primary bulk validation summary; retain individual cohort/stratum estimates.
3. Add the GSE92415 Mayo/age-adjusted model and E-MTAB tissue/drug-adjusted model as sensitivity analyses.
4. Add a quantitative benchmark subsection reporting AUCs, paired bootstrap confidence intervals and gene overlap.
5. Keep the myeloid module as the primary module; describe barrier maturity as a secondary replicated direction; keep IFN/damage and stromal claims cautious.
6. State explicitly that no additional independent single-cell anti-TNF outcome cohort was found.
7. Do not add GSE73661, GSE12251 or GSE14580 as independent cohorts.

## Output map

- `random_effects_meta_analysis.csv`: primary and sensitivity meta-analyses.
- `cross_cohort_signature_effects.csv`: all cohort/stratum effect sizes and AUCs.
- `GSE92415_disease_activity_adjusted_models.csv`: Mayo/age-adjusted models.
- `E-MTAB-7604_tissue_treatment_adjusted_models.csv`: tissue/drug-adjusted models.
- `GSE92415_paired_auc_comparisons.csv` and `E-MTAB-7604_paired_auc_comparisons.csv`: paired benchmark tests.
- `module_vs_published_signature_overlap.csv`: overlap/Jaccard/hypergeometric audit.
- `candidate_dataset_audit.csv` and `covariate_availability_audit.csv`: inclusion and confounder-control audit trails.
- `scripts/run_enhancement_gate.py`: reproducible analysis script.

No manuscript or submission-package file was modified during this gate.
