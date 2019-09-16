# Table of Contents
1. [VETA - Variant prEdiction Tools evAluation](#veta)
2. [Installation](#installation)
3. [Scores annotation](#scores-annotation)
4. [Running on reference datasets](#reference)
5. [Running on unlabelled VCFs](#unlabelled)

<a name="veta"></a>
## VETA - Variant prEdiction Tools evAluation
VETA is a framework that analyses the performance of several variant prediction methods at different levels. It can:
  * 1) Analyse the tools prediction scores on an unseen VCF taking into account the reference thresholds described in literature.
  * 2) Given labeled variants (benign and pathogenic), evaluate tools performance producing several metrics and plots.
  * 3) Inspect reliability of reference thresholds using Clinvar database
  * 4) Apply machine learning to combine scores and improve predictions

<a name="installation"></a>
## Intallation
VETA requires python3 and conda/pip is recommended to install all the dependencies. You can

~~~~
conda create -n veta python=3.6
conda activate veta
git clone https://github.com/PedroBarbosa/VETA.git
conda install --yes --file requirements.conda
pip install -r requirements.pip
~~~~

Alternatively, there is a docker image with all dependencies pre-installed:
~~~~
docker pull pbarbosa/veta:latest
docker run pbarbosa/prediction_tools_evaluation:latest python /tools/VETA/src/veta.py --help
~~~~
To use VETA on your own data (local directory with benign and pathogenic variants OR a single VCF), please map your local data onto the container by seting up a bindmount volume
~~~~
docker run -it -v /local/dir/with/data/:/media pbarbosa/veta:latest
root@e1d195af4858:/tools# veta --dataset /media [Other options..]
~~~~



<a name="scores-annotation"></a>
## Scores annotation

Currently, hg19 version is used, as some tools do not provide scores for the latest genome build. We apply a combination of VEP, 
dbNSFP and vcfanno to annotate variants with the scores. 

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

<a name="reference"></a>
## Running on reference datasets

By default, Clinvar database is used as the ground truth dataset. The simplest command will benchmark prediction tools and produce plots and statistics about their performance. 

`python veta.py -o out_dir`

By default, Clinvar database is filtered so that only variants with 3 stars (and pathogenic, likely_pathhogenic, benign, likely_benign interpretations) are used. If one wants to evaluate Clinvar variants with a different review status, you can set the `--clinvarStars` argument to be less / more stringent in the selection of the ground truth data. Available values for this argument can checked with the `--listClinvarLevels` flag. To run VETA over the clinvar data with 2 stars:

`python veta.py -o out_dir --clinvarStars 2s`

By default, threshold analysis is turned off. To enable it, set the `--thresholdAnalysis` flag. It will use the Clinvar 3 stars with likely assignments (3s_l) as the ground truth data regardless of the argument set in `--clinvarStars`.

Different datasets (other than clinvar) can be used for analysis. Swissvar, Humsavar and Unifun are available. To use a different one, set the `--dataset` argument.

Machine learning analysis can be enbabled with the `--machineLearning` flag. Tools scores will be combined and different standard classifiers are tested.


<a name="unlabelled"></a>
## Running on unlabelled VCFs
