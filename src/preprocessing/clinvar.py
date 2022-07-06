import logging
import os
import sys
import pandas as pd
from preprocessing.location import *

def remove_clinvar_useless(df: pd.DataFrame):
    """
    Removes conflicting, missing and
    uncertain interpretations

    :param pd.DataFrame df: Input df
    :return pd.DataFrame: Df with variants
    following the target significance values
    """
    return df.loc[df['CLNSIG'].isin(['Pathogenic', 'Benign',
                                     'Likely_pathogenic', 'Likely_benign',
                                     'Pathogenic/Likely_pathogenic', 'Benign/Likely_benign'])]


def get_clinvar_cached(fname: str):
    """
    Reads cached processed clinvar file if
    it exists. Else, processes original VCF

    :param str fname: Path to the input File
    :return pd.DataFrame: Processed clinvar df
    """
    tsv_file = fname if fname.endswith('tsv') else fname + '.tsv'

    #If clinvar was processed before
    if os.path.exists(tsv_file):
        logging.info('Cached clinvar file found. VETA will pick from that.' 
                     'If you want to run veta with different configurations, please delete the {} file and run again.'.format(tsv_file))
        return pd.read_csv(tsv_file, low_memory=False)

    return

def filter_by_condition(df: pd.DataFrame, ids:list):
    """
    Filter Clinvar database by a subset of OMIM/MedGen/MONDO IDs
    
    :param list ids: IDs from different databases.
    
    :return pd.DataFrame: Returns the variants in which at least one 
    ID in one database was found
    """

    if all(x is None for x in ids):
        return df
    
    def _get_matched_ids(df, db_field, db_ids):
        """
        Returns Clinvar IDs that match the given condition ID
        """
        re = "[0-9]+" if db_field == "MONDO:" else "[0-9A-Za-z]+"
       
        _ids = df["CLNDISDB"].str.extractall("({}{})".format(db_field, re)).iloc[:, 0].apply(lambda x: x.replace('{}'.format(db_field), '').rstrip())
        match_idx = _ids[_ids.isin(db_ids)].index.get_level_values(0)
        _df = df.loc[match_idx, :]

        return _df.id.to_list()

    ids_map = {0: 'OMIM:', 1: 'MedGen:', 2: 'MONDO:'}
    out_ids = []

    for i, db_ids in enumerate(ids):
        if db_ids is None:
            continue
        
        db_field = ids_map[i]
        out_ids.extend(_get_matched_ids(df, db_field, db_ids))
 
    df = df[df.id.isin(out_ids)]
    if df.shape[0] < 1:
        raise ValueError('No variants left after filtering by {} condition ID(s).'.format(ids))
    else:
        logging.info('Number of variants after filtering by condition ID(s) ({}): {}'.format([x for x in ids if x is not None], df.shape[0]))
   
    return df 

def filter_clinvar_sure(df: pd.DataFrame):
    """
    Removes variants with likely assignments
    """
    return df[df['CLNSIG'].isin(['Pathogenic', 'Benign'])].copy(deep=True)


def filter_clinvar_1_star(df: pd.DataFrame):
    """
    Removes variants with less than 1 star
    """
    # One submitter provided an interpretation with
    # assertion criteria (criteria provided, single submitter)
    # or multiple submitters provided assertion criteria but there
    # are conflicting interpretations. In such cases
    # the independent values are enumerated for clinical
    # significance (criteria provided, conflicting interpretations)
    to_remove = ['no_assertion_criteria_provided',
                 'no_assertion_provided',
                 'no_interpretation_for_the_single_variant']

    return df[~df['CLNREVSTAT'].isin(to_remove)].copy()


def filter_clinvar_2_stars(df: pd.DataFrame):
    """
    Removes variants with less than 2 stars
    """
    return df[(df["CLNREVSTAT"].str.contains('criteria_provided,_multiple_submitters,_no_conflicts')) |
              (df["CLNREVSTAT"].str.contains('reviewed_by_expert_panel')) |
              (df["CLNREVSTAT"].str.contains('practice_guideline'))].copy()


def filter_clinvar_3_stars(df: pd.DataFrame):
    """
    Removes variants with less than 3 stars
    """

    _df = df[df.CLNREVSTAT.str.contains('reviewed_by_expert_panel')].copy()
    return _df


def filter_clinvar_4_stars(df: pd.DataFrame):
    """
    Removes variants with less than 4 stars
    """
    return df[(df["CLNREVSTAT"].str.contains('practice_guideline'))].copy()

