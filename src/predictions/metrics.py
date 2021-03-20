import logging
import sys
from typing import Union, List

import numpy as np
import pandas as pd

from src.preprocessing.utils_tools import ratio

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict
from src.preprocessing import utils_tools
from sklearn.metrics import auc, average_precision_score, roc_auc_score, roc_curve, precision_recall_curve


def generate_statistics(df: pd.DataFrame,
                        statistics: defaultdict,
                        filter_location: str,
                        tool: str,
                        is_single_label: bool = False,
                        f_beta: Union[List, int] = None):
    """
    Compute all stats for a single tool

    :param pd.DataFrame df: Input df with predictions
        compared to the reference threshold
    :param defaultdict statistics: Dict with outputs stats
        to be updated
    :param str filter_location: Location filter
    :param str tool: Tool to compute stats
    :param bool is_single_label: Whether df refers
        to an analysis in the inspect mode with a
        single label. Default: `False`.
    :param Union[List, int] f_beta: Calculate F beta score for the given
    beta value(s). Default: `None`
    :return defaultdict: Updated dict
    """
    # s_df = df[~df[tool + '_prediction'].isnull()]
    total = df.shape[0]
    s_df = df[~df[tool].isnull()]
    statistics['filter'].append(filter_location)
    statistics['tool'].append(tool)

    if is_single_label is False:
        tp = np.sum(s_df['label'].eq(s_df[tool + '_prediction']) & s_df['label'])
        tn = np.sum(s_df['label'].eq(s_df[tool + '_prediction']) & ~s_df['label'])

        fp = np.sum(s_df['label'].ne(s_df[tool + '_prediction']) & ~s_df['label'])
        fn = np.sum(s_df['label'].ne(s_df[tool + '_prediction']) & s_df['label'])

        mp = np.sum(df[tool + '_prediction'].isnull() & df['label'])
        mn = np.sum(df[tool + '_prediction'].isnull() & ~df['label'])
        nan = np.sum(df[tool + '_prediction'].isnull())
        correct = tp + tn

        precision = utils_tools.ratio(tp, tp + fp)
        recall = utils_tools.ratio(tp, tp + fn)
        coverage = utils_tools.ratio(tp + tn + fp + fn, total)
        accuracy = utils_tools.ratio(correct, (total - nan))

        statistics['precision'].append(precision)
        statistics['specificity'].append(utils_tools.ratio(tn, tn + fp))
        statistics['sensitivity'].append(recall)
        statistics['tp'].append(tp)
        statistics['fp'].append(fp)
        statistics['tn'].append(tn)
        statistics['fn'].append(fn)
        statistics['mp'].append(mp)
        statistics['mn'].append(mn)

        statistics['scored_p'].append(np.sum(s_df['label']))
        statistics['scored_n'].append(np.sum(~s_df['label']))
        statistics['total_p'].append(np.sum(df['label']))
        statistics['total_n'].append(np.sum(~df['label']))

        try:
            statistics['F1'].append(round(2 * (precision * recall) /
                                          (precision + recall), 2))

            statistics['weighted_F1'].append(round((2 * (precision * recall) /
                                                    (precision + recall)) * coverage, 2))

        except ZeroDivisionError:
            statistics['F1'].append(0)
            statistics['weighted_F1'].append(0)

        if f_beta is not None:
            if isinstance(f_beta, int):
                f_beta = [f_beta]

            for _b in f_beta:
                try:
                    statistics['Fbeta_' + str(_b)].append(round((1.0 + _b ** 2) * (precision * recall) /
                                                                ((precision * _b ** 2) + recall), 3))

                    statistics['weighted_Fbeta_' + str(_b)].append(
                        round(coverage * (1.0 + _b ** 2) * (precision * recall) /
                              ((precision * _b ** 2) + recall), 3))
                except ZeroDivisionError:
                    statistics['Fbeta_' + str(_b)].append(0)
                    statistics['weighted_Fbeta_' + str(_b)].append(0)

    else:
        correct = np.sum(s_df['label'].eq(s_df[tool]))
        nan = np.sum(df[tool].isnull())
        accuracy = ratio(correct, (total - nan))
        coverage = ratio((total - nan), total)

    weighted_accuracy = round(accuracy * coverage, 2)
    statistics['total'].append(total)
    statistics['correct'].append(correct)
    statistics['nan'].append(nan)
    statistics['fraction_nan'].append(ratio(nan, total))
    statistics['coverage'].append(coverage)
    statistics['accuracy'].append(accuracy)
    statistics['weighted_accuracy'].append(weighted_accuracy)

    return statistics


def do_roc_analysis(data: pd.DataFrame,
                    name: str):
    """
    Perform ROC analysis a given tool

    :param pd.Dataframe data: Scores for a given tool
    with corresponding labels
    :param str name: Tool name
    """
    max_thr = data[name].max()  # + (df_tool[tool].max()) - df_tool[tool].min()) * 0.001
    min_thr = data[name].min()  # - (df_tool[tool].max()) - df_tool[tool].min()) * 0.001

    if pd.isnull(max_thr):
        logging.info("{} has no predicted variants".format(name))
        return

    elif max_thr == min_thr:
        logging.info("{} has the the same min/max score".format(name))
        return

    fpr, tpr, thresh_roc = roc_curve(data.label, data[name])
    precision, recall, thresh_pr = precision_recall_curve(data.label, data[name])

    label = name + "(n=" + str(data.shape[0]) + ","

    if len(thresh_roc) > 400:
        idx = np.round(np.linspace(0, len(thresh_roc) - 1, 200)).astype(int)
        thresh_roc = thresh_roc[idx]
        tpr = tpr[idx]
        fpr = fpr[idx]

    precision = precision[:-1]
    recall = recall[:-1]
    if len(thresh_pr) > 400:
        idx = np.round(np.linspace(0, len(thresh_pr) - 1, 200)).astype(int)
        thresh_pr = thresh_pr[idx]
        precision = precision[idx]
        recall = recall[idx]

    roc_curve_data = [label, list(thresh_roc), list(tpr), list(fpr)]
    pr_curve_data = [label, list(thresh_pr), list(recall), list(precision)]

    return roc_curve_data, pr_curve_data, roc_auc_score(data.label, data[name]), average_precision_score(data.label, data[name])
