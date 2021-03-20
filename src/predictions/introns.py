import logging
import sys

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from collections import defaultdict
from typing import List
import pandas as pd
from src.plots.plots_intronic_analysis import *
from src.plots.plots_benchmark_mode import plot_allele_frequency, plot_unscored, plot_metrics
from src.predictions import metrics
from src.predictions.apply import apply_tool_predictions


def do_intron_analysis(df: pd.DataFrame, thresholds: List, metric: str,
                       out_dir: str,
                       af_column: str,
                       min_variants: int = 20):
    """
    Perform analysis of intronic variants
    in a bin-based manner

    :param pd.DataFrame df: Input df
    :param List thresholds: List with
        the reference thresholds
    :param str metric: Metric to rank the
        tools in the metrics plot
    :param str out_dir: Output directory
    :param str af_column: Allele frequency column
    :param int min_variants: Minimum number
    of variants required to do ROC analysis.
    Default: `20`
    :return:
    """

    logging.info("-------------------------")
    logging.info("Intronic analysis started")
    logging.info("-------------------------")
    assert "intron_bin" in df.columns, "Intronic bins not in the data. " \
                                       "Probably a first run of Clinvar " \
                                       "was performed without --intronic-bins " \
                                       "args. To fix, just remove the 'tsv' in " \
                                       "the input directory and try again."

    out_dir = os.path.join(out_dir, "intron_analysis")
    os.makedirs(out_dir)

    # intronic variants
    df_i = df[~df['intron_bin'].isnull()].copy()
    df_i = apply_tool_predictions(df_i, thresholds)

    try:
        plot_general_bin_info(df_i, out_dir, af_column)
    except ValueError:
        logging.info("Problem plotting info about intronic bins. Skipping.")

    roc_per_bin = defaultdict(list)
    na = {}
    for _bin, _filter_func in filter_intronic_bins:
        logging.info("Looking at {} bin.".format(_bin))
        _df_i_bin = _filter_func(df_i).copy()
        _df_i_bin['count_class'] = _df_i_bin.groupby('outcome')['outcome'].transform('size')

        if _df_i_bin.shape[0] < min_variants:
            logging.info("Not enough variants for ROC analysis "
                         "at {} bin ({})".format(_bin, _df_i_bin.shape[0]))
            continue

        n_pos = _df_i_bin["label"].sum()
        n_neg = (~ _df_i_bin["label"]).sum()
        logging.info("Pathogenic: {}".format(n_pos))
        logging.info("Benign: {}".format(n_neg))

        list_df_metrics_per_tool, roc_metrics_per_tool, pr_metrics_per_tool, general_metrics_per_tool = [], [], [], []
        statistics = defaultdict(list)
        stats_df = ""
        for tool, direction, threshold, *args in thresholds:

            try:
                _no_null_df = _df_i_bin.loc[pd.notnull(_df_i_bin[tool + "_prediction"]), ].copy()
            except KeyError:
                na[tool] = 1

            statistics = metrics.generate_statistics(_df_i_bin, statistics, _bin, tool)
            stats_df = pd.DataFrame(statistics)

            na[tool] = stats_df.loc[stats_df['tool'] == tool, 'fraction_nan'].iloc[0]
            df_tool = _no_null_df[~_no_null_df[tool + "_prediction"].isnull()]

            ####################
            ### ROC analysis ###
            ####################
            if df_tool.shape[0] > min_variants and min(n_pos, n_neg) > 10:

                roc_curve, pr_curve, roc_auc, ap_score = metrics.do_roc_analysis(df_tool[[tool, 'label']],
                                                                                 tool)
                roc_metrics_per_tool.append([tool, float(na[tool]), roc_curve[0], roc_curve[1],
                                             roc_curve[2], roc_curve[3], roc_auc])
                pr_metrics_per_tool.append([tool, float(na[tool]), pr_curve[0], pr_curve[1],
                                            pr_curve[2], pr_curve[3], ap_score])

            else:
                logging.info("ROC analysis will be skipped for {} at {} bin. Not enough variants were predicted "
                             "(minimum accepted: {}; predicted: {}) or the number of counts in the minority class "
                             "is less than 10 ({} variants).".format(tool, _bin, min_variants,
                                                                     df_tool.shape[0], min(n_pos, n_neg)))

            f1_score = stats_df.loc[stats_df['tool'] == tool, 'F1'].iloc[0]
            weighted_f1_score = stats_df.loc[stats_df['tool'] == tool, 'weighted_F1'].iloc[0]
            if not _bin in ["all_intronic", "all_except_0-2", "all_except_0-10", None]:

                roc_per_bin[tool].append([_bin,
                                          f1_score,
                                          weighted_f1_score,
                                          na[tool]])

        stats_df.drop(['filter'], axis=1).to_csv(os.path.join(out_dir, "statistics_{}.csv".format(_bin)),
                                                 sep="\t",
                                                 index=False)

        af_plot = os.path.join(out_dir, "AF_{}".format(_bin))
        unscored_plot = os.path.join(out_dir, "unscored_fraction_{}".format(_bin))
        metrics_plot = os.path.join(out_dir, "tools_metrics_{}".format(_bin))

        plot_allele_frequency(_df_i_bin, af_plot, af_column)
        plot_unscored(stats_df, unscored_plot)
        plot_metrics(stats_df, metrics_plot, metric)

        if roc_metrics_per_tool:
            plot_curve(roc_metrics_per_tool,
                       os.path.join(out_dir, "{}_ROC".format(_bin)),
                       (n_pos, n_neg))

        if pr_metrics_per_tool:
            plot_curve(pr_metrics_per_tool,
                       os.path.join(out_dir, "{}_ROC_pr".format(_bin)),
                       (n_pos, n_neg),
                       is_roc=False)

    df_roc_bins = pd.DataFrame([[k] + i for k, v in roc_per_bin.items() for i in v],
                               columns=["tool", "bin", "F1", "weighted_F1", "fraction_nan"])

    plot_metrics_by_bin(df_roc_bins, os.path.join(out_dir, "per_bin"))
