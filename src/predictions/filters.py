from typing import List
from collections import defaultdict
import sys
import logging
from preprocessing.location import CONSEQUENCES

def update_thresholds(config_dict: defaultdict):
    """
    Updates the list of reference
    thresholds based in the tools
    config file provided and parsed.
    Returns the available methods that exist
    by default in VETA, as well as custom
    methods provided by the user.

    :return List: Complete list of reference
    threshold
    """
    _default_list = [
        ('GERP', '>', 4.4, 'Conservation'),
        ('phyloP', '>', 1.6, 'Conservation'),
        ('SiPhy', '>', 12.17, 'Conservation'),
        ('phastCons', '>', 0.99, 'Conservation'),
        ('CDTS', '<', 10, 'Conservation'),
        ('fitCons', '>', 0.4, 'Whole_genome'),
        ('LINSIGHT', '>', 0.056, 'Whole_genome'), # S-CAP paper

        ('Sift', '<', 0.001, 'Protein'),
        ('Polyphen2HVAR', '>', 0.978, 'Protein'),
        ('Polyphen2HDIV', '>', 0.978, 'Protein'),
        ('MutationAssessor', '>', 1.935, 'Protein'),
        ('MutationTaster2', '>', 0.5, 'Protein'),
        ('FATHMM', '<', -4.14, 'Protein'),
        ('Provean', '<', -2.5, 'Protein'),
        ('LRT', '<', 0.01, 'Protein'),
        ('Mutpred', '>', 0.5, 'Protein'),
        ('VEST4', '>', 0.764, 'Protein'),
        ('CAROL', '>', 0.98, 'Protein'),
        ('Condel', '>', 0.468, 'Protein'),
        ('REVEL', '>', 0.644, 'Protein'),
        ('MetaLR', '>', 0.5, 'Protein'),
        ('MetaSVM', '>', 0.5, 'Protein'),
        ('M-CAP', '>', 0.025, 'Protein'),
        ('MVP', '>', 0.7, 'Protein'),
        ('MTR', '<', 0.5, 'Protein'), # Did not find in the paper, assumed 0.5
        ('MPC', '>', 1.36, 'Protein'),
        ('MISTIC', '>', 0.5, 'Protein'),
        ('CardioBoost', '>', 0.9, 'Protein'),
        ('CardioVAI', '>', 2, 'Protein'),
        ('PrimateAI', '>', 0.790, 'Protein'),
        ('VARITY', '>', 0.75, 'Protein'), # VARITY_R model
        ('ClinPred', '>', 0.5, 'Protein'),
        ('MutScore', '>', 0.5, 'Protein'), # Did not find in the paper, assumed 0.5
        ('MutFormer', '>', 0.5, 'Protein'), # Did not find in the paper, assumed 0.5

        ('cVEP', '>', 0.5, 'Protein'), # 0.5 as artificial threshold
        ('EVE_class20', '>', 0.5, 'Protein'), # 0.5 as artificial threshold
        ('EVE_class50', '>', 0.5, 'Protein'), # 0.5 as artificial threshold
        ('EVE_class90', '>', 0.5, 'Protein'), # 0.5 as artificial threshold
        ('EVE', '>', 0.5, 'Protein'),

        ('CAPICE', '>', 0.02, 'Whole_genome'),
        ('CADD_v1.5', '>', 15, 'Whole_genome'),
        ('CADD-Splice', '>', 15, 'Whole_genome'),

        ('DANN', '>', 0.9, 'Whole_genome'),
        ('GWAVA', '>', 0.5, 'Whole_genome'),
        ('FATHMM-MKL', '>', 0.5, 'Whole_genome'),
        ('Eigen', '>', 4.86, 'Whole_genome'),  # S-CAP paper
        ('ReMM', '>', 0.984, 'Whole_genome'),
        ('FunSeq2', '>', 1.5, 'Whole_genome'),

        # HAL model scores alt5 PSI. kipoi veff scores the DIFF
        # between the ALT and REF allele, thus I set 5 as the
        # minimum threshold to account for PSI changes caused
        # by variant
        ('HAL', '>', 5, 'Splicing'),
        ('MMSplice', '>', 2, 'Splicing'),
        # ('mmsplice-deltaLogitPSI', '>', 2, 'Splicing'),
        # ('mmsplice-pathogenicity', '>', 0.5, 'Splicing'),
        # ('mmsplice-efficiency', '>', 1, 'Splicing'),
        ('kipoiSplice4', '>', 0.5, 'Splicing'),
        ('kipoiSplice4_cons', '>', 0.5, 'Splicing'),

        # S-CAP model has region specific thresholds.
        # To generalize, if score is greater than its
        # specific threshold, I sum 1 to the score
        # (see process_scap_scores method)
        ('ConSpliceML', '>', 0.5, 'Splicing'),
        ('S-CAP', '>', 1, 'Splicing'),
        # TraP, same as S-CAP
        ('TraP', '>', 1, 'Splicing'),
        ('SpliceAI', '>', 0.2, 'Splicing'),
        # SPIDEX: dpsi zscore; dpsi_maxtissue >=5
        # is used in the paper.
        ('SPIDEX', '>', 2, 'Splicing'),
        # Both ada_score and rf_score must be higher
        # than 0.6. Check dbscSNV_merge method on how
        # I dealtwith exceptions
        ('dbscSNV', '>', 0.6, 'Splicing'),
        # MaxEntScan - naive estimation of prediction score.
        # It refers to those variants where the difference in
        # the maximum entropy estimation between the REF and ALT
        # for splice site detection is higher than 3, be it for
        # the gain (positive entropy value) or for the loss
        # (negative entropy value).
        ('MaxEntScan', '>', 3, 'Splicing'),
        ('SQUIRLS', '>', 0.074, 'Splicing'), #CI-SpliceAI paper
        ('IntSplice2', '>', 0.5, 'Splicing'),
        ('CI-SpliceAI', '>', 0.190, 'Splicing'),
        ('Pangolin', '>', 0.2, 'Splicing'),
        ('SPiP', '>', 0.452, 'Splicing'),
        ('MLCsplice', '>', 0.5, 'Splicing'),
        ('AbSplice-DNA', '>', 0.01, 'Splicing'),
        ('LaBranchoR', '>', 0.1, 'Splicing')
    ]

    DEFAULT_TOOLS = [t[0] for t in _default_list]
    _updated_list = []
    for _tool, info in config_dict.items():

        # If user provided custom thresholds
        # for tools available in VETA by
        # default
        if _tool in DEFAULT_TOOLS and len(info) > 1:
            logging.info("Warn: Threshold information about "
                         "{} tool already exists by default "
                         "in VETA. This threshold will be "
                         "updated to the one provided in this run:"
                         "\'{}\'.".format(_tool, info))

            updated_info = tuple([_tool] + info[1:])
            _updated_list.append(updated_info)

        # Use default setting for a
        # known tool
        elif _tool in DEFAULT_TOOLS:
            _default_info = [v for v in _default_list if v[0] == _tool][0]
            _updated_list.append(_default_info)

        # If a custom method
        else:
            if len(info) <= 4:
                custom = [_tool] + info[1:]
            # If custom function, remove it
            elif len(info) == 5:
                custom = [_tool] + info[1:-1]
                
            _updated_list.append(tuple(custom))

    return _updated_list

def _extract_possible_filters(aggregate: bool):
    """
    Extracts possible location/class filters
    to perform an analysis, given the aggregate
    option provided
    
    :return list: List with lambda function to filter df by each loc
    """
    idx = 0 if aggregate else 1
    all_locs = list(set([x[idx] for x in CONSEQUENCES.values()]))
    all_locs.append('deep_intronic')
    all_locs.insert(0, 'all')
    return all_locs

def subset_variants_by_type(types_to_analyse: List = None):
    """
    Subset filters for variant types
    :param List types_to_analyse: List of types
    to filter. Default: `None`
    :return List: Variant types that will be used
    """
    if types_to_analyse is None:
        return filters_var_type
    else:
        return [var_types for var_types in filters_var_type if var_types[0] in types_to_analyse]


def subset_toolset_by_scope(threshold_list: List,
                            scopes: List = None):

    """Returns subset of tools belonging to the
    scope defined, or if `is_intronic` is true,
    removes tools belonging to `Protein scope,
    regardless of the `scores_to_analys` value

    :param List threshold_list: List of reference
        thresholds defined taking into account the
        default thresholds in VETA as well as the custom
        ones provided by the user

    :param List scopes: Tool scopes to analyse.
        Default: `None`. Use all scopes
    """

    _to_analyse = ['Conservation', 'Protein', 'Whole_genome', 'Splicing']

    if scopes is not None:
        _to_analyse = list(set(_to_analyse).intersection(scopes))

    out = []
    for tool in threshold_list:
        try:
            if tool[3] in _to_analyse:
                out.append(tool)
        except IndexError:
            raise IndexError('Scope is missing in the config file for {} tool'.format(tool))
    return out


filters_var_type = [
    ('all_types', lambda x: x),
    ('snps', lambda x: x[x['type'].str.match('snp')]),
    ('indels', lambda x: x.query('type == "indel" & subtype == "ins" | type == "indel" & subtype == "del"')),
    ('insertions', lambda x: x.query('type == "indel" & subtype == "ins"')),
    ('deletions', lambda x: x.query('type == "indel" & subtype == "del"')),
    ('mnps', lambda x: x.query('type == "indel" & subtype == "mnp"'))
]

filter_intronic_bins = [
    ('all_intronic', lambda x: x[~x['intron_bin'].isnull()]),
    ('all_except_1-2', lambda x: x[~x['intron_bin'].str.match('1-2')]),
    ('all_except_1-10', lambda x: x[(~x['intron_bin'].str.match('1-2')) &
                                    (~x['intron_bin'].str.match('3-10')) & 
                                    (~x['intron_bin'].str.match('1-10'))]),
    ('1-2', lambda x: x[x['intron_bin'].str.match('1-2')]),
    ('3-10', lambda x: x[x['intron_bin'].str.match('3-10')]),
    ('1-10', lambda x: x[x['intron_bin'].str.match('1-10')]),
    ('11-40', lambda x: x[x['intron_bin'].str.match('11-40')]),
    ('41-200', lambda x: x[x['intron_bin'].str.match('41-200')]),
    ('201-500', lambda x: x[x['intron_bin'].str.match('201-500')]),
    ('501-1000', lambda x: x[x['intron_bin'].str.match('501-1000')]),
    ('1000+', lambda x: x[x['intron_bin'].str.match('1000\+')])
]