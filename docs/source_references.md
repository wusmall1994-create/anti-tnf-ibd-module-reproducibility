# Source references and dataset citations

## Primary single-cell resource

1. **GSE282122 / TAURUS atlas** - GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282122  
   Notes: Public on November 19, 2024; 216 gut biopsy libraries; title: *A longitudinal single-cell atlas of anti-tumour necrosis factor treatment in inflammatory bowel disease*.

2. **TAURUS Nature Immunology article** - PubMed: https://pubmed.ncbi.nlm.nih.gov/39438660/ ; Nature: https://www.nature.com/articles/s41590-024-01994-8  
   Notes: reports ~1 million single-cell transcriptomes from 216 gut biopsies and 41 subjects following adalimumab therapy.

3. **TAURUS Zenodo release** - https://zenodo.org/records/14007626  
   Notes: annotated H5AD lineage objects are hosted here. The 2026-07-26 lineage-resolved update of this package uses these annotated H5AD objects (see `docs/taurus_zenodo_acquisition_2026-07-26.md`); the earlier GEO-matrix workflow does not require them.

4. **TAURUS public code** - https://github.com/DendrouLab/TAURUS_paper

## Bulk anti-TNF validation cohorts

5. **GSE16879** - GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16879  
   Notes: mucosal expression profiling in IBD before and after first infliximab treatment; 133 samples on GPL570.

6. **GSE14580** - GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14580 ; PubMed: https://pubmed.ncbi.nlm.nih.gov/19700435/  
   Notes: UC mucosal gene signatures for infliximab response; appears to overlap with the UC subset of GSE16879 in this analysis and should not be counted independently.

7. **GSE23597** - GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE23597 ; related PubMed article: https://pubmed.ncbi.nlm.nih.gov/21448149/  
   Notes: colonic biopsy expression data from infliximab/placebo-treated UC patients; 113 samples from 48 patients at baseline, week 8 and week 30.

8. **GSE12251** - GEO accession page: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12251  
   Notes: predictive response signature to infliximab in UC; downloaded/available as candidate validation but not yet deeply analysed in this project.

## Additional biological context included in the final reference list

9. **Anti-TNF resistance cellular module context** - Martin JC et al. Cell 2019: https://www.sciencedirect.com/science/article/pii/S0092867419308967

10. **OSM/OSMR anti-TNF response context** - West NR et al. Nat Med 2017: https://pmc.ncbi.nlm.nih.gov/articles/PMC5420447/

11. **TREM1 anti-TNF response context** - Verstockt B et al. EBioMedicine 2019: https://www.sciencedirect.com/science/article/pii/S2352396419300325

12. **Cell-centred meta-analysis of anti-TNF non-response** - Gaujoux R et al. Gut 2019: https://gut.bmj.com/content/68/4/604

13. **Systematic review of predictive biomarkers** - Stevens TW et al. Aliment Pharmacol Ther 2018: https://www.ovid.com/journals/alpt/fulltext/10.1111/apt.15033~systematic-review-predictive-biomarkers-of-therapeutic

## Citation caution

For the final manuscript, verify bibliographic fields from PubMed or journal pages and export to RIS/BibTeX before submission. GEO pages should be cited for data accessions; original articles should be cited for study design and biological interpretation.
