import logging
import numpy as np
import pandas as pd
from typing import Union, List
from collections import defaultdict
from preprocessing.utils_tools import ratio
from sklearn.metrics import auc, average_precision_score, roc_auc_score, roc_curve, precision_recall_curve, matthews_corrcoef


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

        precision = ratio(tp, tp + fp)
        recall = ratio(tp, tp + fn)
        coverage = ratio(tp + tn + fp + fn, total)
        accuracy = ratio(correct, (total - nan))
        
        #mcc = ((tp * tn) - (fp * fn)) / (sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
        mcc = round(matthews_corrcoef(s_df.label.astype(int), s_df[tool + '_prediction'].astype(int)), 3)
        normalized_mcc = (mcc + 1) / 2
        weighted_norm_mcc = round(normalized_mcc * coverage, 3)
        
        statistics['precision'].append(precision)
        statistics['specificity'].append(ratio(tn, tn + fp))
        statistics['sensitivity'].append(recall)
        statistics['mcc'].append(mcc)
        statistics['norm_mcc'].append(normalized_mcc)
        statistics['weighted_norm_mcc'].append(weighted_norm_mcc)
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
                                          (precision + recall), 3))

            statistics['weighted_F1'].append(round((2 * (precision * recall) /
                                                    (precision + recall)) * coverage, 3))

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

    weighted_accuracy = round(accuracy * coverage, 3) 
    statistics['total'].append(total)
    statistics['correct'].append(correct)
    statistics['nan'].append(nan)
    statistics['fraction_nan'].append(ratio(nan, total))
    statistics['coverage'].append(coverage)
    statistics['accuracy'].append(accuracy)
    statistics['weighted_accuracy'].append(weighted_accuracy)

    return statistics


def do_roc_analysis(_data: pd.DataFrame,
                    name: str,
                    higher_is_better: bool = True):
    """
    Perform ROC analysis a given tool

    :param pd.Dataframe data: Scores for a given tool
    with corresponding labels
    :param str name: Tool name
    :param bool higher_is_better: If pathogenicity is predicted when prediction is higher than reference threshold
    """
    
    # ROC analysis is not performed for 
    # S-CAP (it uses multiple reference thresholds) and 
    # cVEP (predictions are discrete categories that VETA converts to fixed floats)

    if name.lower() in ['scap', 's-cap', 'cvep']:
        return

    data = _data.copy()
    max_thr = data[name].max()  # + (df_tool[tool].max()) - df_tool[tool].min()) * 0.001
    min_thr = data[name].min()  # - (df_tool[tool].max()) - df_tool[tool].min()) * 0.001

    if pd.isnull(max_thr):
        logging.info("{} has no predicted variants".format(name))
        return

    elif max_thr == min_thr:
        logging.info("{} has the the same min/max score".format(name))
        return

    if higher_is_better is False:
        # Reverse scores, so that higher is better
        data['y_pred'] = data[name] * (-1)

    else:
        data['y_pred'] = data[name]
    
    data = data[['label', 'y_pred']].sort_values('y_pred')
    data['rank'] = np.arange(len(data))
    data['rank_score'] = data['rank'] / len(data)

    fpr, tpr, thresh_roc = roc_curve(data.label, data.rank_score)
    precision, recall, thresh_pr = precision_recall_curve(data.label, data.rank_score)

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

    return roc_curve_data, pr_curve_data, auc(fpr, tpr), auc(recall, precision)
