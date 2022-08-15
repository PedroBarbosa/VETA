import logging
import pandas as pd
import multiprocessing
from functools import partial
from itertools import repeat
from collections import defaultdict
from tqdm import tqdm
from typing import List
from preprocessing.osutils import ensure_folder_exists
from plots.plots_intronic_analysis import *
from plots.plots_benchmark_mode import (
    plot_allele_frequency,
    plot_unscored,
    plot_metrics,
    plot_curve,
)
from predictions import metrics
from predictions.apply import apply_tool_predictions


def _parallel_roc(thresholds: list, _df_i_bin: pd.DataFrame, **kwargs):

    roc_metrics, pr_metrics = [], []
    statistics = defaultdict(list)
    na = {}
    tool, direction = thresholds[0], thresholds[1]

    try:
        _no_null_df = _df_i_bin.loc[
            pd.notnull(_df_i_bin[tool + "_prediction"]),
        ].copy()
    except KeyError:
        na[tool] = 1

    statistics = metrics.generate_statistics(
        _df_i_bin, statistics, kwargs["_bin"], tool
    )
    stats_df = pd.DataFrame(statistics)

    na[tool] = stats_df.loc[stats_df["tool"] == tool, "fraction_nan"].iloc[0]
    df_tool = _no_null_df[~_no_null_df[tool + "_prediction"].isnull()]

    ####################
    ### ROC analysis ###
    ####################
    if kwargs["skip_roc"] is False:
        n_pos_pred = df_tool[df_tool[tool + "_prediction"]].shape[0]
        n_neg_pred = df_tool[df_tool[tool + "_prediction"] == False].shape[0]

        if float(na[tool]) > 0.5:
            logging.info(
                "ROC analysis will be skipped for {} at {} bin. "
                "More than 50% of missing data ({})".format(
                    tool, kwargs["_bin"], round(float(na[tool]), 2)
                )
            )
      
        elif min(n_pos_pred, n_neg_pred) < 10:
            logging.info(
                "ROC analysis will be skipped for {} at {} bin. "
                "No minimum number of predictions on each class (10) found (N={})".format(
                    tool, kwargs["_bin"], min(n_pos_pred, n_neg_pred)
                )
            )

        else:
            try:
                roc_curve, pr_curve, roc_auc, ap_score = metrics.do_roc_analysis(
                    df_tool[[tool, "label"]], tool, higher_is_better=direction == ">"
                )
                roc_metrics = [
                    tool,
                    float(na[tool]),
                    roc_curve[0],
                    roc_curve[1],
                    roc_curve[2],
                    roc_curve[3],
                    roc_auc,
                ]
                pr_metrics = [
                    tool,
                    float(na[tool]),
                    pr_curve[0],
                    pr_curve[1],
                    pr_curve[2],
                    pr_curve[3],
                    ap_score,
                ]

            # E.g. S-CAP
            except TypeError:
                pass

    f1_score = stats_df.loc[stats_df["tool"] == tool, "F1"].iloc[0]
    weighted_f1_score = stats_df.loc[stats_df["tool"] == tool, "weighted_F1" ].iloc[0]

    return stats_df, roc_metrics, pr_metrics, [tool, f1_score, weighted_f1_score, na[tool]]


def do_intron_analysis(
    df: pd.DataFrame,
    thresholds: List,
    metric: str,
    aggregate_classes: str,
    out_dir: str,
    af_column: str,
    min_variants: int = 20,
):
    """
    Perform analysis of intronic variants
    in a bin-based manner

    :param pd.DataFrame df: Input df
    :param List thresholds: List with
        the reference thresholds
    :param bool aggregate_classes:
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
    assert "intron_bin" in df.columns, (
        "Intronic bins are not present in the data."
        "Probably a first run of Clinvar "
        "was performed without --do_intronic_analysis "
        "args. To fix, just remove the 'tsv' in "
        "the input directory and try again."
    )

    out_dir = os.path.join(out_dir, "intron_analysis")
    out_tsv = os.path.join(out_dir, "results_tsv")
    out_af = os.path.join(out_dir, "allele_frequency")
    out_fixed_thresh = os.path.join(out_dir, "performance_at_fixed_thresh")
    out_roc = os.path.join(out_dir, "roc_analysis")
    out_all_bin_agg = os.path.join(out_dir, "all_bin_together")
    [
        ensure_folder_exists(path)
        for path in [
            out_dir,
            out_tsv,
            out_af,
            out_fixed_thresh,
            out_roc,
            out_all_bin_agg,
        ]
    ]

    # intronic variants
    df_i = df[~df["intron_bin"].isnull()].copy(deep=True)
    df_i = apply_tool_predictions(df_i, thresholds)

    metrics_per_bin = defaultdict(list)
    loc_in_protein_coding = [
        "splice_site",
        "splice_region",
        "intronic",
        "deep_intronic",
    ]
    pc_is_done = False
    bin_to_exclude = ["1-2", "3-10"] if aggregate_classes else ["1-10"]
    locations = sorted(df_i.location.unique().tolist())
    locations.insert(0, "all")

    for _class in locations:

        if _class == "all":
            _df_loc = df_i.copy(deep=True)
        elif _class in loc_in_protein_coding:
            if pc_is_done is False:
                _class = "no_UTRs"
                _df_loc = df_i[df_i.location.isin(loc_in_protein_coding)].copy(
                    deep=True
                )
                pc_is_done = True
            else:
                continue
        else:
            _df_loc = df_i[df_i.location == _class].copy(deep=True)

        print("\n")
        logging.info("--------------")
        logging.info("Looking at {} intronic variants".format(_class))
        logging.info("--------------")
        print()
        if _df_loc.empty:
            logging.info("No {} variants found.".format(_class))
            continue

        _df_loc["count_class"] = _df_loc.groupby("outcome")["outcome"].transform("size")

        try:
            plot_general_bin_info(
                _df_loc, out_dir, _class, aggregate_classes, af_column
            )
        except ValueError:
            logging.info("Problem plotting info about intronic bins. Skipping.")

        for _bin, _filter_func in filter_intronic_bins:
            if "all_except" in _bin:
                continue

            if _bin in bin_to_exclude:
                continue
            skip_roc = False
            print()
            logging.info("-------")
            logging.info("BIN: {}.".format(_bin))
            logging.info("-------")
            _df_i_bin = _filter_func(_df_loc).copy(deep=True)

            if _df_i_bin.empty:
                logging.info("No {} variants at {} bin.".format(_class, _bin))
                continue

            n_pos = _df_i_bin["label"].sum()
            n_neg = (~_df_i_bin["label"]).sum()
            logging.info("Pathogenic: {}".format(n_pos))
            logging.info("Benign: {}".format(n_neg))

            if _df_i_bin.shape[0] < min_variants:
                logging.info(
                    "Not enough {} variants for ROC analysis "
                    "at {} bin (N={}, required={})".format(
                        _class, _bin, _df_i_bin.shape[0], min_variants
                    )
                )
                skip_roc = True

            elif min(n_pos, n_neg) < 10:
                logging.info(
                    "Number of variants in the minority class is lower (N={}) "
                    "than minimum required (N=10) for ROC analysis.".format(
                        min(n_pos, n_neg)
                    )
                )
                skip_roc = True

            kwargs = {
                "_class": _class,
                "_bin": _bin,
                "skip_roc": skip_roc,
            }

            with multiprocessing.Pool() as p:

                stats_df, roc_metrics, pr_metrics, _per_bin = zip(*p.map(partial(_parallel_roc, 
                                                                  _df_i_bin=_df_i_bin,
                                                                  **kwargs),
                                                          thresholds))
                stats_df = pd.concat(stats_df)
 
            if not _bin in ["all_intronic", "all_except_1-2", "all_except_1-10", None]:
                for tool_info in _per_bin:
                    metrics_per_bin[tool_info[0]].append([_class, _bin, 
                                                      tool_info[1],
                                                      tool_info[2], 
                                                      tool_info[3]])
   
            stats_out = os.path.join(
                out_tsv, "statistics_{}_{}.tsv".format(_class, _bin)
            )
            stats_df.drop(["filter"], axis=1).to_csv(stats_out, sep="\t", index=False)
            af_plot = os.path.join(out_af, "AF_{}_{}.pdf".format(_class, _bin))
            unscored_plot = os.path.join(
                out_fixed_thresh, "unscored_fraction_{}_{}.pdf".format(_class, _bin)
            )
            metrics_plot = os.path.join(
                out_fixed_thresh, "tools_metrics_{}_{}.pdf".format(_class, _bin)
            )

            plot_allele_frequency(_df_i_bin, af_plot, af_column)
            plot_unscored(stats_df, unscored_plot)
            plot_metrics(stats_df, metrics_plot, metric)

            plot_curve(
                    roc_metrics,
                    os.path.join(out_roc, "{}_{}_ROC.pdf".format(_class, _bin)),
                    (n_pos, n_neg),
                )

            plot_curve(
                    pr_metrics,
                    os.path.join(out_roc, "{}_{}_ROC_pr.pdf".format(_class, _bin)),
                    (n_pos, n_neg),
                    is_roc=False,
                )

        df_metrics_per_bin = pd.DataFrame(
            [[k] + i for k, v in metrics_per_bin.items() for i in v],
            columns=[
                "tool",
                "variant_class",
                "bin",
                "F1",
                "weighted_F1",
                "fraction_nan",
            ],
        )

        df_metrics_per_bin.groupby("variant_class").apply(
            plot_metrics_by_bin,
            os.path.join(out_all_bin_agg, "per_bin"),
            aggregate_classes,
        )
