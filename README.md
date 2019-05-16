# VETA - Variant prEdiction Tools evAluation

VETA is a framework that analyses the performance of several variant prediction methods at different levels. It can:
  * 1) Analyse the tools prediction scores on an unseen VCF taking into account the reference thresholds described in literature.
  * 2) Given labeled variants (benign and pathogenic), evaluate tools performance producing several metrics and plots.
  * 3) Inspect reliability of reference thresholds using Clinvar database
  * 4) Apply machine learning to combine scores and improve predictions

### Scores ###

Currently, hg19 version is used, as some tools do not provide scores for the latest genome build. We apply a combination of VEP, 
dbNSFP and vcfanno to annotate variants with the scores. We expect to release in a near future a simple web application to distance
users from the cumbersome process of running the tools and/or annotating VCFs with the scores. 

---
List of tools available (more will be continuosly added)
---
| Category                        | Tool              | Threshold            | Incorporates_other_scores | Obtained_from |
|---------------------------------|-------------------|----------------------|---------------------------|---------------|
| Conservation scores             | GERP              | > 4.4                | No                        | vcfanno       |
|                                 | phyloP            | > 1.6                | No                        | vcfanno       |
|                                 | SiPhy             | > 12.7               | No                        | vcfanno       |
|                                 | phastCons         | > 0.8                | No                        | vcfanno       |
| Effects on protein function     | SIFT              | < 0.05               | No                        | dbNSFP v4.02  |
|                                 | PolyPhen-2 HDIV   | > 0.5                | No                        | dbNSFP v4.02  |
|                                 | PolyPhen-2 HVAR   | > 0.5                | No                        | dbNSFP v4.02  |
|                                 | LRT               | < 0.01               | No                        | dbNSFP v4.02  |
|                                 | Mutation Assessor | > 1.935              | No                        | dbNSFP v4.02  |
|                                 | FATHMM            | < -1.5               | No                        | dbNSFP v4.02  |
|                                 | PROVEAN           | < -2.5               | No                        | dbNSFP v4.02  |
|                                 | MutationTaster2   | > 0.5                | No                        | dbNSFP v4.02  |
|                                 | MutPred           | > 0.5                | Yes                       | dbNSFP v4.02  |
|                                 | Condel            | >  0.98              | Yes                       | VEP plugin    |
|                                 | CAROL             | > 0.468              | Yes                       | VEP plugin    |
|                                 | M-CAP             | > 0.025              | Yes                       | dbNSFP v4.02  |
|                                 | REVEL             | > 0.5                | Yes                       | dbNSFP v4.02  |
|                                 | VEST4             | > 0.67               | Yes                       | dbNSFP v4.02  |
|                                 | MetaSVM           | > 0.5                | Yes                       | dbNSFP v4.02  |
|                                 | MetaLR            | > 0.5                | Yes                       | dbNSFP v4.02  |
| Functional genomics annotations | FATHMM-MKL        | > 0.5                | Yes                       | vcfanno       |
|                                 | GWAVA             | > 0.4                | Yes                       | vcfanno       |
|                                 | Eigen             | > 1                  | Yes                       | vcfanno       |
|                                 | ReMM              | > 0.984              | Yes                       | custom_script |
| Fitness measurement             | CADD v1.4         | > 15                 | Yes                       | VEP plugin    |
|                                 | DANN              | > 0.9                | Yes                       | custom_script |
|                                 | fitCons           | > 0.4                | No                        | vcfanno       |
|                                 | LINSIGHT          | > 0.4                | Yes                       | vcfanno       |
| Splicing                        | MaxEntScan        | \|MaxEntScan_diff\|> 1 | No                        | VEP plugin    |
|                                 | dbscSNV           | > 0.6                | Yes                       | VEP plugin    |
|                                 | SPIDEX            | \|dpsi_zscore\| > 2    | No                        | vcfanno       |
|                                 | traP              | > 0.92               | Yes                       | vcfanno       |
|                                 | SpliceAI          | > 0.2                | No                        | vcfanno       |
