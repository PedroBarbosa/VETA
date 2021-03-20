from typing import List
from collections import defaultdict
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


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
        ('fitCons', '>', 0.4, 'Functional'),
        ('LINSIGHT', '>', 0.4, 'Functional'),

        ('Sift', '<', 0.05, 'Protein'),
        ('Polyphen2HVAR', '>', 0.5, 'Protein'),
        ('Polyphen2HDIV', '>', 0.5, 'Protein'),
        ('MutationAssessor', '>', 1.935, 'Protein'),
        ('MutationTaster', '>', 0.5, 'Protein'),
        ('FATHMM', '<', -1.5, 'Protein'),
        ('Provean', '<', -2.5, 'Protein'),
        ('LRT', '<', 0.01, 'Protein'),
        ('Mutpred', '>', 0.5, 'Protein'),
        ('VEST4', '>', 0.67, 'Protein'),
        ('CAROL', '>', 0.98, 'Protein'),
        ('Condel', '>', 0.468, 'Protein'),
        ('REVEL', '>', 0.5, 'Protein'),
        ('MetaLR', '>', 0.5, 'Protein'),
        ('MetaSVM', '>', 0.5, 'Protein'),
        ('M-CAP', '>', 0.025, 'Protein'),
        ('MVP', '>', 0.7, 'Protein'),
        ('CardioBoost', '>', 0.9, 'Protein'),
        ('PrimateAI', '>', 0.8, 'Protein'),

        ('CADD', '>', 15, 'Functional'),
        ('DANN', '>', 0.9, 'Functional'),
        ('GWAVA', '>', 0.5, 'Functional'),
        ('FATHMM-MKL', '>', 0.5, 'Functional'),
        ('Eigen', '>', 1, 'Functional'),  # M-CAP paper
        ('ReMM', '>', 0.984, 'Functional'),
        ('FunSeq2', '>', 1.5, 'Functional'),

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
        ('MaxEntScan', '>', 3, 'Splicing')
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
                         "in VETA. These thresholds will be "
                         "updated to the one provided by the user:"
                         "\'{}\'. Take caution with such setup "
                         "because the default thresholds in VETA "
                         "should be fine as they were carefully "
                         "extracted from literature".format(_tool, info))

            updated_info = tuple([_tool] + info[1:])
            _updated_list.append(updated_info)

        # Use default setting for a
        # known tool
        elif _tool in DEFAULT_TOOLS:
            _default_info = [v for v in _default_list if v[0] == _tool][0]
            _updated_list.append(_default_info)

        # If a custom method
        else:
            custom = [_tool] + info[1:]
            _updated_list.append(custom)

    return _updated_list

filters_location = [
    ('all', lambda x: x),
    ('coding', lambda x: x[x['location'].str.match('coding')]),
    ('splice_site', lambda x: x[x['location'].str.match('splice_site')]),
    ('splice_region', lambda x: x[x['location'].str.match('splice_region')]),
    ('5primeUTR', lambda x: x[x['location'].str.match('5primeUTR')]),
    ('3primeUTR', lambda x: x[x['location'].str.match('3primeUTR')]),
    ('deep_intronic', lambda x: x[x['location'].str.match('deep_intronic')]),
    ('noncodingRNA', lambda x: x[x['location'].str.match('noncodingRNA')]),
    ('mithocondrial', lambda x: x[x['location'].str.match('mithocondrial')]),
    ('unknown', lambda x: x[x['location'].str.match('unknown')])
]

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
    ('all_except_0-2', lambda x: x[~x['intron_bin'].str.match('0-2')]),
    ('all_except_0-10', lambda x: x[~x['intron_bin'].str.match('0-10')]),
    ('0-2', lambda x: x[x['intron_bin'].str.match('0-2')]),
    ('3-10', lambda x: x[x['intron_bin'].str.match('3-10')]),
    ('11-30', lambda x: x[x['intron_bin'].str.match('11-30')]),
    ('31-100', lambda x: x[x['intron_bin'].str.match('31-100')]),
    ('101-200', lambda x: x[x['intron_bin'].str.match('101-200')]),
    ('201-500', lambda x: x[x['intron_bin'].str.match('201-500')]),
    ('500+', lambda x: x[x['intron_bin'].str.match('500\+')])
]


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
                            scopes: List = None,
                            is_intronic: bool = False):

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
    :param bool is_intronic: Whether analysis
        is targeted for intronic variants.
        Default: `False`. If `True`, 'Protein'
        scope will not be used.
    """

    _to_analyse = ['Conservation', 'Protein', 'Functional', 'Splicing']

    if is_intronic:
        _to_analyse.remove('Protein')

    if scopes is not None:
        _to_analyse = list(set(_to_analyse).intersection(scopes))

    try:
        return [tool for tool in threshold_list if tool[3] in _to_analyse]
    except IndexError:
        raise IndexError('Scope is missing in the config file for some custom provided tool')
