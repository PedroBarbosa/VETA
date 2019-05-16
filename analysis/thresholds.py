threshold_list = [
    ('GERP', '>', 4.4,'#00441b','o'),
    ('phyloP', '>', 0.5,'#00441b',','), # TODO
    ('SiPhy', '>', 12.16,'#00441b','p'),
    ('phastCons', '>', 0.8,'#00441b','v'), # TODO
    ('MutationAssessor', '>', 1.935,'#00441b','^'),
    ('fitCons', '>', 0.4, '#00441b','h'),
    ('LINSIGHT', '>', 0.4, '#00441b', 'h'),

    ('Sift', '<', 0.05,'#084594','o'),
    ('Polyphen2HVAR', '>', 0.5,'#084594',','),
    ('Polyphen2HDIV', '>', 0.5,'#084594','p'),
    ('MutationTaster', '>', 0.5,'#084594','v'),
    ('FATHMM-MKL', '>', 0.5,'#084594','^'),
    ('FATHMM', '<', -1.5,'#084594','<'),
    ('Provean', '<', -2.5,'#084594','>'),
    # Identification of deleterious mutations within three human genomes
    ('LRT', '<', 0.001,'#084594','h'),
    # Sift-based. Is it 1 or 2? Check out mutpred splice
    ('Mutpred', '>', 0.5,'#084594','8'), #previously 0.611
    ('VEST4','>', 0.67,'#084594','D'),

    # Probabililty of a gene being loss-of-function intolerant (pLI)
    #('ExACpLI', '>', 0.9,'#084594','h'),

    ('CADD', '>', 15,'#7f0000','o'),
    ('DANN', '>', 0.9,'#7f0000',','),
    ('MetaLR', '>', 0.5,'#7f0000','p'),
    ('MetaSVM', '>', 0.5,'#7f0000','v'),
    ('Revel', '>', 0.5,'#7f0000','^'),
    ('GWAVA', '>', 0.5,'#7f0000','<'),
    ('Eigen', '>', 1,'#7f0000','>'), #M-CAP paper
    #('Eigen-PC', '>', 1),  #Not found
    ('ReMM', '>', 0.984,'#7f0000','h'),
    ('CAROL','>', 0.98,'#7f0000','8'),
    ('Condel','>', 0.468,'#7f0000','D'),
    ('M-CAP','>', 0.025,'#7f0000','H'),

    #Splicing related scores
    ('traP', '>', 0.92, '#4a1486','o'), #probably damaging the transcript >= 0.93
    ('SpliceAI', '>', 0.2, '#4a1486','o'),
    ('SPIDEX', '>', 2,'#4a1486','o'), #dpsi zscore; dpsi_maxtissue >=5 is used in the paper.
    ('dbscSNV', '>', 0.6,'#4a1486',','), #both ada_score and rf_score must be higher than 0.6. Check dbscSNV_merge method on how I dealt with exceptions
    ('MaxEntScan','>',1,'#4a1486','p') #naive estimation of prediction score. It refers to those variants where the difference in the maximum entropy estimation between the REF and ALT
    #  for splice site detection is higher than 1, be it for the gain (positive entropy value) or for the loss (negative entropy value). May not be ideal
    #('SpliceRegion','>',0.9) #splice region VEP plugin just detects and report if a variant is located within a splice region, thus this fixed value
]
