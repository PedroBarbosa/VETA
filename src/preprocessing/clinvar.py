import sys
import os
import logging
import numpy as np
import pandas as pd
import hgvs.parser
from .location import get_location, get_location_from_consequence
from .vcf_cleaning import vcf_cleaning
from .vcf_processing import process_vcf_scores
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


def tuple2float(x):
    if not x:
        return np.nan
    return x[1] if x[1] != '.' or x[1] != '' or x[1] else np.nan


def fix_columns(df):
    new_col_names = ['hg19.chr', 'hg19.pos', 'ref', 'alt', 'id', 'type', 'subtype', 'rsID', 'HGVSc', 'Gene', 'Consequence']
    for column in df:
        if isinstance(df[column].iloc[0], (tuple,)):
            new_col_names.append(df[column].iloc[0][0])
            df[column] = df[column].map(tuple2float)

    rename_dict = {i: j for i, j in zip(list(df), new_col_names)}
    df.rename(columns=rename_dict, inplace=True)
    return df


def remove_clinvar_useless(df):
    """ Removes conflicting, missing and uncertain interpretations """
    return df.loc[df['CLNSIG'].isin(['Pathogenic', 'Benign', 'Likely_pathogenic', 'Likely_benign'])]


def get_clinvar_cached(fname, loc):
    tsv_file = fname + ".tsv"
    if os.path.exists(tsv_file):
        return pd.read_csv(tsv_file)
    else:
        print("Clinvar 'tsv' file not found. Dataset will be transformed")
        df = get_clinvar(fname, loc)
        df.to_csv(tsv_file)
        print("Done")
        return df


def get_clinvar(fname, loc):

    scores = process_vcf_scores(fname, is_clinvar=True)
    df = pd.DataFrame.from_dict(scores, orient='index')
    df = fix_columns(df)
    df = df.dropna(subset=["CLNSIG","CLNREVSTAT"])
    df = remove_clinvar_useless(df)

    if loc == "HGVSc":
        hp = hgvs.parser.Parser()
        df['location'] = df['HGVSc'].apply(get_location,hp=hp)
    else:
        df['location'] = df['Consequence'].apply(get_location_from_consequence)

    df['class'] = df['CLNSIG'].str.contains('Pathogenic')    
    df = vcf_cleaning(df)
    return df.reset_index()


def filter_clinvar_sure(df):
    return df[df['CLNSIG'].isin(['Pathogenic', 'Benign'])].copy()


def filter_clinvar_1_star(df):
    #	One submitter provided an interpretation with assertion criteria (criteria provided, single submitter)
    #  or multiple submitters provided assertion criteria but there are conflicting interpretations in which case
    # the independent values are enumerated for clinical significance (criteria provided, conflicting interpretations)
    #print(df["clinvar.rcv.review_status"].value_counts())
    df1 = df[(~df["CLNREVSTAT"].str.contains('no_assertion_criteria_provided'))
              | (~df["CLNREVSTAT"].str.contains('no_assertion_provided'))
              | (~df["CLNREVSTAT"].str.contains('no_interpretation_for_the_single_variant'))]
    return df1.copy()


def filter_clinvar_2_stars(df):
    df1 = df[(df["CLNREVSTAT"].str.contains('criteria_provided,_multiple_submitters,_no_conflicts')) |
              (df["CLNREVSTAT"].str.contains('reviewed_by_expert_panel')) |
              (df["CLNREVSTAT"].str.contains('practice_guideline'))]
    return df1.copy()


def filter_clinvar_3_stars(df):
    df1 = df[(df["CLNREVSTAT"].str.contains('reviewed_by_expert_panel'))]
    return df1.copy()    


def filter_clinvar_4_stars(df):
    df1 = df[(df["CLNREVSTAT"].str.contains('practice_guideline'))]
    return df1.copy()