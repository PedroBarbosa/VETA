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

        logging.info("Pathogenic: {}".format(_df_i_bin["label"].sum()))
        logging.info("Benign: {}".format((~ _df_i_bin["label"]).sum()))

        list_df_metrics_per_tool = []
        statistics = defaultdict(list)
        stats_df = ""
        for tool, direction, threshold, *args in thresholds:

            try:
                _no_null_df = _df_i_bin.loc[pd.notnull(_df_i_bin[tool + "_prediction"]),].copy()
            except KeyError:
                na[tool] = 1
                continue

            statistics = metrics.generate_statistics(_df_i_bin, statistics, _bin, tool)
            stats_df = pd.DataFrame(statistics)

            na[tool] = stats_df.loc[stats_df['tool'] == tool, 'fraction_nan'].iloc[0]
            df_tool = _no_null_df[~_no_null_df[tool + "_prediction"].isnull()]
            if df_tool.shape[0] > min_variants:

                # Do ROC
                tool_metrics, roc_auc, pr_auc = metrics.do_roc_analysis(df_tool[[tool, 'label']],
                                                                        tool,
                                                                        direction)

                f1_score = [stats_df.loc[stats_df['tool'] == tool, 'F1'].iloc[0]] * len(tool_metrics)
                weighted_f1_score = [stats_df.loc[stats_df['tool'] == tool, 'weighted_F1'].iloc[0]] * len(tool_metrics)
                nan = [na[tool]] * len(tool_metrics)

                tool_metrics = [x + [y] + [z] + [f] + [wf] + [na] for x, y, z, f, wf, na in zip(tool_metrics, roc_auc,
                                                                                                pr_auc,
                                                                                                f1_score,
                                                                                                weighted_f1_score,
                                                                                                nan)]
                list_df_metrics_per_tool.append(pd.DataFrame(tool_metrics,
                                                             columns=["tool", "threshold",
                                                                      "precision", "recall",
                                                                      "FPR", "ROC-auc", "PR-auc",
                                                                      "F1", "weighted_F1",
                                                                      "fraction_nan"]))

                if not _bin in ["all_intronic", "all_except_0-10", None]:
                    roc_per_bin[tool].append([_bin,
                                              roc_auc[0],
                                              pr_auc[0],
                                              f1_score[0],
                                              weighted_f1_score[0],
                                              na[tool]]
                                             )
            else:
                logging.info("Not enough variants predicted by {} "
                             "for ROC analysis at {} bin".format(tool,
                                                                 _bin))
                continue

        stats_df.drop(['filter'], axis=1).to_csv(os.path.join(out_dir, "statistics_{}.csv".format(_bin)),
                                                 sep="\t",
                                                 index=False)

        af_plot = os.path.join(out_dir, "AF_{}".format(_bin))
        unscored_plot = os.path.join(out_dir, "unscored_fraction_{}".format(_bin))
        metrics_plot = os.path.join(out_dir, "tools_metrics_{}".format(_bin))

        plot_allele_frequency(_df_i_bin, af_plot, af_column)
        plot_unscored(stats_df, unscored_plot)
        plot_metrics(stats_df, metrics_plot, metric)

        final_df_metrics = pd.concat(list_df_metrics_per_tool)

        try:
            plot_ROCs(final_df_metrics,
                      os.path.join(out_dir, "{}_ROC".format(_bin)),
                      _df_i_bin.loc[_df_i_bin['outcome'] == "Pathogenic", 'count_class'].iloc[0]
                      )
        except IndexError:
            logging.info("No positive class instances at {} bin, skipping ROC.".format(_bin))

    df_roc_bins = pd.DataFrame([[k] + i for k, v in roc_per_bin.items() for i in v],
                               columns=["tool", "bin", "auROC", "prROC",
                                        "F1", "weighted_F1",  "fraction_nan"])

    plot_metrics_by_bin(df_roc_bins, os.path.join(out_dir, "per_bin_evol"))