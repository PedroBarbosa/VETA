threshold_list_complete = [
    ('GERP', '>', 4.4, 'Conservation', '#00441b','o'),
    ('phyloP', '>', 0.5, 'Conservation', '#00441b',','), # TODO
    ('SiPhy', '>', 12.16, 'Conservation', '#00441b','p'),
    ('phastCons', '>', 0.8, 'Conservation', '#00441b','v'), # TODO
    ('fitCons', '>', 0.4, 'Functional', '#00441b','h'),
    ('LINSIGHT', '>', 0.4, 'Functional', '#00441b', 'h'),

    ('Sift', '<', 0.05, 'Protein', '#084594','o'),
    ('Polyphen2HVAR', '>', 0.5, 'Protein', '#084594',','),
    ('Polyphen2HDIV', '>', 0.5, 'Protein', '#084594','p'),
    ('MutationAssessor', '>', 1.935, 'Protein','#00441b','^'),
    ('MutationTaster', '>', 0.5, 'Protein', '#084594','v'),
    ('FATHMM', '<', -1.5, 'Protein', '#084594','<'),
    ('Provean', '<', -2.5, 'Protein', '#084594','>'),
    # Identification of deleterious mutations within three human genomes
    ('LRT', '<', 0.001, 'Protein', '#084594','h'),
    # Sift-based. Is it 1 or 2? Check out mutpred splice
    ('Mutpred', '>', 0.5, 'Protein', '#084594','8'), #previously 0.611
    ('VEST4','>', 0.67, 'Protein', '#084594','D'),
    ('CAROL', '>', 0.98, 'Protein', '#7f0000', '8'),
    ('Condel', '>', 0.468, 'Protein', '#7f0000', 'D'),
    ('Revel', '>', 0.5, 'Protein', '#7f0000', '^'),
    ('MetaLR', '>', 0.5, 'Protein', '#7f0000', 'p'),
    ('MetaSVM', '>', 0.5, 'Protein', '#7f0000', 'v'),
    ('M-CAP', '>', 0.025, 'Protein', '#7f0000', 'H'),
    # Probabililty of a gene being loss-of-function intolerant (pLI)
    #('ExACpLI', '>', 0.9,'#084594','h'),

    ('CADD', '>', 15, 'Functional', '#7f0000','o'),
    ('DANN', '>', 0.9, 'Functional', '#7f0000',','),
    ('GWAVA', '>', 0.5, 'Functional', '#7f0000','<'),
    ('FATHMM-MKL', '>', 0.5, 'Functional', '#084594','^'),
    ('Eigen', '>', 1, 'Functional', '#7f0000','>'), #M-CAP paper
    #('Eigen-PC', '>', 1),  #Not found
    ('ReMM', '>', 0.984, 'Functional', '#7f0000','h'),
    ('funseq2', '>', 1.5, 'Functional', '#7f0000','h'), #paper


    #Splicing related scores
    ('HAL', '>1', 5, 'Splicing', '#4a1486', 'o'), #model scores alt5 PSI. kipoi veff scores the DIFF between the ALT and REF allele, thus we set 5 as the minimum threshold to account for PSI changes caused by variant
    ('MMSplice', '>', 1, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-deltaLogitPSI', '>', 1, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-pathogenicity', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-efficiency', '>', 1, 'Splicing', '#4a1486', 'o'),
    ('kipoiSplice-4', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    ('kipoiSplice-4-cons', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    ('S-CAP', '>', 1, 'Splicing', '#4a1486', 'o'), #model has region specific thresholds. To generalize, if score is greater than its specific threshold, we sum up 1 to the score (see process_scap_scores method)
    ('traP', '>', 0.92, 'Splicing', '#4a1486','o'), #probably damaging the transcript >= 0.93
    ('SpliceAI', '>', 0.2, 'Splicing', '#4a1486','o'),
    ('SPIDEX', '>', 2, 'Splicing', '#4a1486','o'), #dpsi zscore; dpsi_maxtissue >=5 is used in the paper.
    ('dbscSNV', '>', 0.6, 'Splicing', '#4a1486',','), #both ada_score and rf_score must be higher than 0.6. Check dbscSNV_merge method on how I dealt with exceptions
    ('MaxEntScan','>',1, 'Splicing', '#4a1486','p') #naive estimation of prediction score. It refers to those variants where the difference in the maximum entropy estimation between the REF and ALT
    #  for splice site detection is higher than 1, be it for the gain (positive entropy value) or for the loss (negative entropy value). May not be ideal
    #('SpliceRegion','>',0.9) #splice region VEP plugin just detects and report if a variant is located within a splice region, thus this fixed value
]


def subset_toolset_by_scope(thresholds, scope_to_analyse=None):

    if scope_to_analyse:
        return [tool for tool in thresholds if any(True for t in tool[3].split("|") if t in scope_to_analyse)]
    else:
        return thresholds
