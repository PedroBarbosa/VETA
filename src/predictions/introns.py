import logging
import pandas as pd
import multiprocessing
from functools import partial
from collections import defaultdict
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
                    tool, kwargs["_bin"], round(float(na[tool]), 3)
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
    weighted_f1_score = stats_df.loc[stats_df["tool"] == tool, "weighted_F1"].iloc[0]

    return (
        stats_df,
        roc_metrics,
        pr_metrics,
        [tool, f1_score, weighted_f1_score, na[tool]],
    )


class IntronicAnalysis(object):
    """
    Base class for intronic analysis
    """

    def __init__(
        self,
        df: pd.DataFrame,
        split_splice_sites: bool,
        thresholds: List,
        metric: str,
        aggregate_classes: str,
        out_dir: str,
        af_column: str,
    ):
        """
        Perform analysis of intronic variants
        in a bin-based manner

        :param pd.DataFrame df: Input df
        :param bool split_splice_sites: 
        :param List thresholds: List with
            the reference thresholds
        :param bool aggregate_classes:
        :param str metric: Metric to rank the
            tools in the metrics plot
        :param str out_dir: Output directory
        :param str af_column: Allele frequency column
        :return:
        """
        self.split_splice_sites = split_splice_sites
        self.thresholds = thresholds
        self.metric = metric
        self.aggregate_classes = aggregate_classes
        self.out_dir = os.path.join(out_dir, "intron_analysis")
        self.af_column = af_column
        self.metrics_per_bin = defaultdict(list)
        
        assert "intron_bin" in df.columns, (
            "Intronic bins are not present in the data."
            "Probably a first run of Clinvar "
            "was performed without --do_intronic_analysis "
            "args. To fix, just remove the 'tsv' in "
            "the input directory and try again."
        )
        # intronic variants
        self.df = df[~df["intron_bin"].isnull()].copy(deep=True)
        self.df = apply_tool_predictions(self.df, self.thresholds)
        logging.info("-------------------------")
        logging.info("Intronic analysis started")
        logging.info("-------------------------")

        loc_in_protein_coding = [
            "splice_site",
            "splice_region",
            "intronic",
            "deep_intronic",
        ]
        pc_is_done = False
        self.bin_to_exclude = ["1-2", "3-10"] if self.aggregate_classes else ["1-10"]
        locations = sorted(self.df.location.unique().tolist())
        locations.insert(0, "all")

        for _class in locations:
            self._class = _class
            ss = None
            if self._class == "all":
                _df_loc = self.df.copy(deep=True)

                if self.split_splice_sites:
                    out_per_bin_metrics = []
                    for ss in ["donor", "acceptor"]:
                        _df = _df_loc[_df_loc.which_ss == ss].copy()
                        out_per_bin_metrics.append(self.per_bin_computations(_df, specific_ss=ss))

                    plot_metrics_by_bin_split_ss(pd.concat(out_per_bin_metrics),
                                                  os.path.join(self.out_dir, "all_bin_together/per_bin"),
                                                  self.aggregate_classes)
        
                    self.donor_vs_acceptor(
                        _df_loc, os.path.join(self.out_dir, "performance_at_fixed_thresh")
                    )
                    ss=None
                self.per_bin_computations(_df_loc, specific_ss=ss)
                    
            # elif self._class in loc_in_protein_coding:
            #     if pc_is_done is False:
            #         _class = "no_UTRs"
            #         _df_loc = self.df[self.df.location.isin(loc_in_protein_coding)].copy(
            #             deep=True
            #         )
            #         pc_is_done = True
            #     else:
            #         continue
            #     self.per_bin_computations(_df_loc, specific_ss=ss)
            # else:
            #      _df_loc = self.df[self.df.location == _class].copy(deep=True)  
            #     self.per_bin_computations(_df_loc, specific_ss=ss)

    def donor_vs_acceptor(self, df: pd.DataFrame, out_dir: str):

        stats_d, stats_a = defaultdict(list), defaultdict(list)
        out_file = os.path.join(out_dir, "donor_vs_acceptor_{}.pdf".format(self._class))

        do_donor, do_acceptor = True, True
        donor = df[df.which_ss == "donor"]
        acceptor = df[df.which_ss == "acceptor"]

        if all(x.shape[0] < 20 for x in [donor, acceptor]):
            logging.info(
                "Not enough intronic variants (N < 30) assigned as donor or acceptor-associated (< 1000bp from splice site). Skipping."
            )
            return

        elif donor.shape[0] < 20:
            do_donor = False
            logging.info(
                "Not enough donor-associated variants (N = {}). Only acceptor-related variants will be analyzed.".format(
                    donor.shape[0]
                )
            )

        elif acceptor.shape[0] < 20:
            do_acceptor = False
            logging.info(
                "Not enough donor-associated variants (N = {}). Only acceptor-related variants will be analyzed.".format(
                    acceptor.shape[0]
                )
            )

        for tool, _, _, _ in self.thresholds:
            if do_donor:
                stats_d = metrics.generate_statistics(donor, stats_d, "all", tool)

            if do_acceptor:
                stats_a = metrics.generate_statistics(acceptor, stats_a, "all", tool)

        stats_donor = pd.DataFrame(stats_d)
        stats_acceptor = pd.DataFrame(stats_a)

        if all(x is True for x in [do_donor, do_acceptor]):
            n_donors = stats_donor["total"].max()
            n_acceptors = stats_acceptor["total"].max()
            title = "N donors = {}; N acceptors = {}".format(n_donors, n_acceptors)

            perf = pd.merge(
                stats_donor[["tool", self.metric]],
                stats_acceptor[["tool", self.metric]],
                on="tool",
                how="left",
            )
            perf.columns = ["Tool", "Donor", "Acceptor"]
            perf = perf.melt(
                id_vars="Tool",
                var_name="Splice site",
                value_name=self.metric.replace("_", " "),
            )

        elif do_donor:
            n_donors = stats_donor["total"].max()
            title = "N donors = {}".format(n_donors)

            perf = stats_donor[["tool", self.metric]]
            perf.columns = ["Tool", "Donor"]
            perf = perf.melt(
                id_vars="Tool",
                var_name="Splice site",
                value_name=self.metric.replace("_", " "),
            )

        elif do_acceptor:

            n_acceptors = stats_acceptor["total"].max()
            title = "N acceptors = {}".format(n_acceptors)

            perf = stats_acceptor[["tool", self.metric]]
            perf.columns = ["Tool", "Donor"]
            perf = perf.melt(
                id_vars="Tool",
                var_name="Splice site",
                value_name=self.metric.replace("_", " "),
            )

        plot_donor_vs_acceptor(perf, self.metric.replace("_", " "), out_file, title)

    def per_bin_computations(self, 
                             _df_loc, 
                             specific_ss=None, 
                             min_variants=20):
        """
        :param int min_variants: Minimum number
        of variants required to do ROC analysis.
        Default: `20`
        """
        out_tsv = os.path.join(self.out_dir, "results_tsv")
        out_af = os.path.join(self.out_dir, "allele_frequency")
        out_fixed_thresh = os.path.join(self.out_dir, "performance_at_fixed_thresh")
        out_roc = os.path.join(self.out_dir, "roc_analysis")
        out_all_bin_agg = os.path.join(self.out_dir, "all_bin_together")

        [
            ensure_folder_exists(path)
            for path in [
                self.out_dir,
                out_tsv,
                out_af,
                out_fixed_thresh,
                out_roc,
                out_all_bin_agg,
            ]
        ]

        _aux_class = (
            self._class
            if specific_ss is None
            else self._class + "_" + specific_ss + "_related"
        )
        print("\n")
        logging.info("--------------")
        logging.info("Looking at {} intronic variants".format(_aux_class))
        logging.info("--------------")
        print()

        if _df_loc.empty:
            logging.info("No {} variants found.".format(_aux_class))
            return

        _df_loc["count_class"] = _df_loc.groupby("outcome")["outcome"].transform("size")

        try:
            plot_general_bin_info(
                _df_loc,
                self.out_dir,
                _aux_class,
                self.aggregate_classes,
                self.af_column,
            )
        except ValueError:
            logging.info("Problem plotting info about intronic bins. Skipping.")

        for _bin, _filter_func in filter_intronic_bins:
            #if "all_except" in _bin:
            #    continue

            if _bin in self.bin_to_exclude:
                continue
            skip_roc = False
            print()
            logging.info("-------")
            logging.info("BIN: {}.".format(_bin))
            logging.info("-------")
            _df_i_bin = _filter_func(_df_loc).copy(deep=True)

            if _df_i_bin.empty:
                logging.info("No {} variants at {} bin.".format(_aux_class, _bin))
                continue

            n_pos = _df_i_bin["label"].sum()
            n_neg = (~_df_i_bin["label"]).sum()
            logging.info("Pathogenic: {}".format(n_pos))
            logging.info("Benign: {}".format(n_neg))
                
            if _df_i_bin.shape[0] < min_variants:
                logging.info(
                    "Not enough {} variants for ROC analysis "
                    "at {} bin (N={}, required={})".format(
                        _aux_class, _bin, _df_i_bin.shape[0], min_variants
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
                "_bin": _bin,
                "skip_roc": skip_roc,
            }

            with multiprocessing.Pool() as p:

                stats_df, roc_metrics, pr_metrics, _per_bin = zip(
                    *p.map(
                        partial(_parallel_roc, _df_i_bin=_df_i_bin, **kwargs),
                        self.thresholds,
                    )
                )
                stats_df = pd.concat(stats_df)

            roc_m = plot_curve(
                roc_metrics,
                os.path.join(out_roc, "{}_{}_ROC.pdf".format(_aux_class, _bin)),
                (n_pos, n_neg),
            )

            pr_m = plot_curve(
                pr_metrics,
                os.path.join(out_roc, "{}_{}_ROC_pr.pdf".format(_aux_class, _bin)),
                (n_pos, n_neg),
                is_roc=False,
            )

            stats_out = os.path.join(out_tsv, "statistics_{}_{}.tsv".format(_aux_class, _bin))
            af_plot = os.path.join(out_af, "AF_{}_{}.pdf".format(_aux_class, _bin))
            unscored_plot = os.path.join(out_fixed_thresh, "unscored_fraction_{}_{}.pdf".format(_aux_class, _bin))
            metrics_plot = os.path.join(out_fixed_thresh, "tools_metrics_{}_{}.pdf".format(_aux_class, _bin))

            plot_allele_frequency(_df_i_bin, af_plot, self.af_column)
            plot_unscored(stats_df, unscored_plot)
            plot_metrics(stats_df, metrics_plot, self.metric)
            
            if not _bin in ["all_intronic", "all_except_1-2", "all_except_1-10", None]:

                for tool_info in _per_bin:
                    res = [_aux_class, _bin, tool_info[1], tool_info[2], tool_info[3]]

                    if skip_roc:
                        res.extend([None, None])
                    else:
                        res.append(roc_m[tool_info[0]]) if tool_info[0] in roc_m.keys() else res.append(None)
                        res.append(pr_m[tool_info[0]]) if tool_info[0] in pr_m.keys() else res.append(None)
                        
                    self.metrics_per_bin[tool_info[0]].append(res)
            
            if roc_m:
                stats_df = pd.merge(stats_df, pd.DataFrame.from_dict(roc_m, orient='index', columns=['auROC']), how = 'left', left_on = 'tool', right_index = True)
            else:
                stats_df['auROC'] = None
            
            if pr_m:
                stats_df = pd.merge(stats_df, pd.DataFrame.from_dict(pr_m, orient='index', columns=['average_precision']), how = 'left', left_on = 'tool', right_index = True)
            else:
                stats_df['average_precision'] = None

            stats_df.drop(["filter"], axis=1).to_csv(stats_out, sep="\t", index=False)
            
        df_metrics_per_bin = pd.DataFrame(
            [[k] + i for k, v in self.metrics_per_bin.items() for i in v],
            columns=[
                "tool",
                "variant_class",
                "bin",
                "F1",
                "weighted_F1",
                "fraction_nan",
                "roc_auc",
                "ap_score"
            ],
        )

        df_metrics_per_bin.groupby("variant_class").apply(
            plot_metrics_by_bin,
            os.path.join(out_all_bin_agg, "per_bin"),
            self.aggregate_classes,
        )
        
        if specific_ss:
            _df = df_metrics_per_bin[df_metrics_per_bin.variant_class.str.contains('related')]
            return _df
