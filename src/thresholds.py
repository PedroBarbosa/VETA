threshold_list_complete = [
    ('GERP', '>', 4.4, 'Conservation', '#00441b', 'o'),
    ('phyloP', '>', 1.6, 'Conservation', '#00441b', ','),
    ('SiPhy', '>', 12.17, 'Conservation', '#00441b', 'p'),
    ('phastCons', '>', 0.99, 'Conservation', '#00441b', 'v'),
    ('fitCons', '>', 0.4, 'Functional', '#00441b', 'h'),
    ('LINSIGHT', '>', 0.4, 'Functional', '#00441b', 'h'),

    ('Sift', '<', 0.05, 'Protein', '#084594','o'),
    ('Polyphen2HVAR', '>', 0.5, 'Protein', '#084594',','),
    ('Polyphen2HDIV', '>', 0.5, 'Protein', '#084594','p'),
    ('MutationAssessor', '>', 1.935, 'Protein','#00441b','^'),
    ('MutationTaster', '>', 0.5, 'Protein', '#084594','v'),
    ('FATHMM', '<', -1.5, 'Protein', '#084594','<'),
    ('Provean', '<', -2.5, 'Protein', '#084594','>'),
    ('LRT', '<', 0.01, 'Protein', '#084594','h'),

    ('Mutpred', '>', 0.5, 'Protein', '#084594','8'),
    ('VEST4','>', 0.67, 'Protein', '#084594','D'),
    ('CAROL', '>', 0.98, 'Protein', '#7f0000', '8'),
    ('Condel', '>', 0.468, 'Protein', '#7f0000', 'D'),
    ('REVEL', '>', 0.5, 'Protein', '#7f0000', '^'),
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
    #('Eigen', '>', 10, 'Functional', '#7f0000','>'),
    ('ReMM', '>', 0.984, 'Functional', '#7f0000','h'),
    ('FunSeq2', '>', 1.5, 'Functional', '#7f0000','h'),
    #('FunSeq2', '>', 1, 'Functional', '#7f0000','h'),

    #Splicing related scores
    ('HAL', '>', 5, 'Splicing', '#4a1486', 'o'), #model scores alt5 PSI. kipoi veff scores the DIFF between the ALT and REF allele, thus we set 5 as the minimum threshold to account for PSI changes caused by variant
    ('MMSplice', '>', 2, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-deltaLogitPSI', '>', 2, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-pathogenicity', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    #('mmsplice-efficiency', '>', 1, 'Splicing', '#4a1486', 'o'),
    ('kipoiSplice4', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    ('kipoiSplice4_cons', '>', 0.5, 'Splicing', '#4a1486', 'o'),
    ('S-CAP', '>', 1, 'Splicing', '#4a1486', 'o'), #model has region specific thresholds. To generalize, if score is greater than its specific threshold, we sum up 1 to the score (see process_scap_scores method)
    ('TraP', '>', 1, 'Splicing', '#4a1486','o'), #same as scap
    ('SpliceAI', '>', 0.2, 'Splicing', '#4a1486','o'),
    ('SPIDEX', '>', 2, 'Splicing', '#4a1486','o'), #dpsi zscore; dpsi_maxtissue >=5 is used in the paper.
    ('dbscSNV', '>', 0.6, 'Splicing', '#4a1486',','), #both ada_score and rf_score must be higher than 0.6. Check dbscSNV_merge method on how I dealt with exceptions
    ('MaxEntScan', '>', 3, 'Splicing', '#4a1486','p') #naive estimation of prediction score. It refers to those variants where the difference in the maximum entropy estimation between the REF and ALT
    #  for splice site detection is higher than 1, be it for the gain (positive entropy value) or for the loss (negative entropy value). May not be ideal
    #('SpliceRegion','>',0.9) #splice region VEP plugin just detects and report if a variant is located within a splice region, thus this fixed value
]


def subset_toolset_by_scope(thresholds, scope_to_analyse=None, to_intronic=None):
    """Returns subset of tools belonging to the scope defined,
    or if to_intronic is true, removes tools belonging
    to Protein scope, regardless of the scopes_to_analyse value
    """
    if to_intronic:
        return [tool for tool in thresholds if any(True for t in tool[3].split("|") if t != "Protein")]
    elif scope_to_analyse:
        return [tool for tool in thresholds if any(True for t in tool[3].split("|") if t in scope_to_analyse)]
    else:
        return thresholds

