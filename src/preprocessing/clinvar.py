import logging
import os
import sys

import pandas as pd

from src.preprocessing.location import *

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def remove_clinvar_useless(df: pd.DataFrame):
    """
    Removes conflicting, missing and
    uncertain interpretations

    :param pd.DataFrame df: Input df
    :return pd.DataFrame: Df with variants
    following the target significance values
    """
    return df.loc[df['CLNSIG'].isin(['Pathogenic', 'Benign',
                                     'Likely_pathogenic', 'Likely_benign'])]


def get_clinvar_cached(fname: str):
    """
    Reads cached processed clinvar file if
    it exists. Else, processes original VCF

    :param str fname: Path to the input File
    :return pd.DataFrame: Processed clinvar df
    """
    tsv_file = fname if fname.endswith('tsv') else fname + '.tsv'

    # If clinvar was processed before
    if os.path.exists(tsv_file):
        return pd.read_csv(tsv_file, low_memory=False)

    return


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

