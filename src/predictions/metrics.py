import pandas as pd
import numpy as np
import sys
import logging
from src.preprocessing.utils import ratio
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict
from src.preprocessing import utils
from sklearn.metrics import auc
from scipy import integrate


def generate_statistics(df: pd.DataFrame,
                        statistics: defaultdict,
                        filter_location: str,
                        tool: str,
                        is_single_label: bool = False):
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

        precision = utils.ratio(tp, tp + fp)
        recall = utils.ratio(tp, tp + fn)
        coverage = utils.ratio(tp + tn + fp + fn, total)
        accuracy = utils.ratio(correct, (total - nan))

        statistics['precision'].append(precision)
        statistics['specificity'].append(utils.ratio(tn, tn + fp))
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


def do_ROC_analysis(data: pd.DataFrame,
                    name: str,
                    direction: str,
                    n_thresholds: int = 100):
    """
    Perform ROC analysis a given tool

    :param pd.Dataframe data: Scores for a given tool
    with corresponding labels
    :param str name: Tool name
    :param str direction: Direction of pathogenicity
        score
    :param int n_thresholds: Number of different
        thresholds to measure. Default: `100`
    """
    max_thr = data[name].max()  # + (df_tool[tool].max()) - df_tool[tool].min()) * 0.001
    min_thr = data[name].min()  # - (df_tool[tool].max()) - df_tool[tool].min()) * 0.001

    if pd.isnull(max_thr):
        logging.info("{} has no predicted variants".format(name))
        return

    elif max_thr == min_thr:
        logging.info("{} has the the same min/max score".format(name))
        return

    _step = (max_thr - min_thr) / float(n_thresholds)
    threshold_range = np.arange(min_thr, max_thr, _step)

    tool_metrics = []

    for _thresh in threshold_range:

        if direction == ">":
            classification_f = lambda x: x == np.nan and np.nan or x >= _thresh
        else:
            classification_f = lambda x: x == np.nan and np.nan or x <= _thresh

        classification = data[name].map(classification_f)
        tp = np.sum(classification.eq(data['label']) & classification)
        fp = np.sum(~data['label'] & classification)
        fn = np.sum(data['label'] & ~classification)
        tn = np.sum(classification.eq(data['label']) & ~classification)
        ap = tp + fn
        ap_predicted = np.sum(classification)

        sensitivity = ratio(tp, ap)
        precision = ratio(tp, ap_predicted)

        an = tn + fp
        an_predicted = np.sum(~classification)
        specificity = ratio(tn, an)

        tool_metrics.append([name + "(n=" + str(data.shape[0]) + ",",
                             _thresh,
                             precision,
                             sensitivity,
                             1 - specificity])

    prc = [val_at_thresh[2] for val_at_thresh in tool_metrics]
    tpr = [val_at_thresh[3] for val_at_thresh in tool_metrics]
    fpr = [val_at_thresh[4] for val_at_thresh in tool_metrics]

    roc_auc = [auc(sorted(fpr), sorted(tpr))] * len(tool_metrics)
    # roc_auc = [np.trapz(tpr, fpr, dx=_step)] * len(tool_metrics)

    pr_auc = [auc(sorted(tpr), sorted(prc))] * len(tool_metrics)

    # pr_auc = [integrate.simps(prc, dx=_step)] * len(tool_metrics)
    # trapz gives the same are, but with negative values.
    # Need to reverse X
    # pr_auc = [np.trapz(prc, tpr)] * len(tool_metrics)

    return tool_metrics, roc_auc, pr_auc
