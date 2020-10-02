from src.preprocessing.osutils import check_file_exists
from typing import List
from src.predictions.filters import filters_location
from src.base import Base
from src.predictions.apply import apply_tool_predictions
from src.predictions.metrics import generate_statistics
from src.preprocessing.osutils import ensure_folder_exists
from src.plots.plots_inspect_mode import *
from collections import defaultdict
import pandas as pd
import numpy as np
import os
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
TOOLS_CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "map_tools2vcf_annotation.txt")


class PredictionsEval(Base):
    """
    Inspect predictions of an unlabelled VCF
    """

    def __init__(self, vcf: str,
                 out_dir: str,
                 scope_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 genome: str = "hg19",
                 is_intronic: bool = False,
                 best_tools: str = None,
                 n_best_tools: int = None,
                 plot_these_tools: List = None,
                 labels: str = None,
                 tools_config: str = TOOLS_CONFIG):

        """
        ----
        Base args described in Base class
        ----
        :param str best_tools: Restrict analysis to the best set
            of tools obtained from a previous VETA run using a
            reference catalog (e.g. Clinvar). It must refer to the
            file `tools_ranking*.csv` that is written when running
            the aforementioned analysis. Default: `None`, use all
            tools  available considering the `scope_to_predict`
            argument.

        :param int n_best_tools: Number of best tools selected
            from the ranking provided in the `--best_tools` argument.
            Default: `5`.

        :param List plot_these_tools: Plot scores distribution for the
            given tools. Default: `None`.

        :param str labels: If variants stored in VCF represent a list
            of labelled variants of a given type. If so, additional
            metrics will be inferred. Available types: ['benign',
            'pathogenic']. Default: `None`, variants are unlabelled.
        """

        super().__init__(vcf=vcf,
                         out_dir=out_dir,
                         scope_to_predict=scope_to_predict,
                         types_of_variant=types_of_variant,
                         metric=metric,
                         location=location,
                         genome=genome,
                         is_intronic=is_intronic,
                         is_clinvar=False,
                         tools_config=tools_config)

        if best_tools is not None:
            check_file_exists(best_tools)

        if plot_these_tools is not None:
            assert set(plot_these_tools).issubset(self.available_tools), "Invalid tool(s) name(s) in the " \
                                                                         "'--plot_these_tools' argument. " \
                                                                         "List of available tools:\n{}".\
                format(self.available_tools)

        if "f1" in self.metric and self.labels is not None:
            raise ValueError("F1-based metrics are not available in the "
                             "inspect mode, since all variants refer to a "
                             "single class type (--label set to {})".format(self.labels))

        self.best_tools = check_file_exists(best_tools)
        self.n_best_tools = n_best_tools
        self.plot_these_tools = plot_these_tools
        self.labels = labels
        self.inspect_predictions()

    def inspect_predictions(self):
        """
        Inspect predictions for each
        variant type and produces output
        files and plots

        :return:
        """

        df_with_pred = apply_tool_predictions(self.df, self.thresholds).set_index('id')

        for var_type, _func in self.variant_types:

            outdir = os.path.join(self.out_dir, var_type)

            ensure_folder_exists(outdir)
            _df_by_type = _func(df_with_pred).copy()

            logging.info("Looking at {} ({} variants)".format(var_type,
                                                              _df_by_type.shape[0]))
            if _df_by_type.shape[0] == 0:
                logging.warning("WARN: There are no {} in the variant set. "
                                "Skipping this analysis.".format(var_type))
                continue

            # df = apply_tool_predictions(df, threshold_list)
            _df_just_pred = pd.concat([_df_by_type[[col for col in _df_by_type.columns if '_prediction' in col]],
                                       _df_by_type["HGVSc"], _df_by_type["location"].to_frame()],
                                      axis=1)

            _df_just_pred = _df_just_pred.rename(columns={col: col.split('_prediction')[0]
                                                          for col in _df_just_pred.columns})
            _df_just_pred.to_csv(os.path.join(outdir, "predictions.tsv"), sep="\t")

            try:
                ratios_df = _df_just_pred.drop(["location", "HGVSc"], axis=1)
                ratios_df = ratios_df.apply(lambda x: x.value_counts(True, dropna=False),
                                            axis=1).fillna(0).sort_values([True], ascending=False)
                ratios_df.rename(columns={False: 'is_benign',
                                          True: 'is_pathogenic',
                                          np.nan: "unpredictable"},
                                 inplace=True)

                plot_area(ratios_df, outdir)
                ratios_df["unpredictable"] *= 100

                plot_heatmap(ratios_df, outdir)
                _top_predicted_pathogenic = ratios_df[ratios_df.is_pathogenic > 0.5]
                if _top_predicted_pathogenic.shape[0] > 0:
                    plot_heatmap(_top_predicted_pathogenic, outdir, display_annot=True)

            except KeyError:
                raise KeyError("No tool has predictions for the given variant type ({}), "
                               "analysis is going to be skipped.".format(var_type))

            # if VCF refers to variants of a
            # given label, compute additional
            # metrics
            if self.labels is not None:
                _df_just_pred['label'] = False if self.labels == "Benign" else True
                self.generate_performance_with_label(_df_just_pred, filters=filters_location,
                                                     thresholds=self.thresholds,
                                                     metric=self.metric,
                                                     outdir=outdir)

            if self.best_tools is not None:
                tools = pd.read_csv(self.best_tools, sep="\t")['tool'].head(self.n_best_tools).tolist()
                tools.append("location")
                _df_top_tools = _df_just_pred[tools]
                class_dict = {True: 1, False: -1}
                plot_heatmap_toptools(_df_top_tools.replace(class_dict),
                                      filters=filters_location,
                                      outdir=outdir)

            if self.plot_these_tools is not None:
                for tool in self.plot_these_tools:
                    # uses original df that contains
                    # scores, not just binary predictions
                    plot_tool_score_distribution(_df_by_type,
                                                 tool=tool,
                                                 thresholds=self.thresholds,
                                                 outdir=outdir)

        if self.is_intronic:
            logging.info("For now, intronic analysis is not available "
                         "in the inspect mode. Skipping it.")

    def generate_performance_with_label(self,
                                        df_pred: pd.DataFrame,
                                        filters: List,
                                        thresholds: List,
                                        metric: str,
                                        outdir: str):
        """
        Evaluate tools performance considering that the
        intput VCF refers to a list of variants with a
        single and known label (e.g. benign or pathogenic)

        :param pd.DataFrame df_pred: Df with predictions
            for each variant
        :param List filters: Location filters to employ
            so that variants at each location are processed
            independently
        :param List thresholds: List of tools with the
            reference thresholds
        :param str metric: Metric to evaluate predictions
        :param str outdir: Output directory
        :return:
        """
        for filter_name, _func in filters:
            outfile = os.path.join(outdir, "tools_ranking_{}.csv").format(filter_name)
            statistics = defaultdict(list)
            df = _func(df_pred).copy()
            if df.shape[0] < 10:
                logging.warning("WARN: Input VCF has not a minimum number of {} "
                                "variants (10) to evaluate tools performance "
                                "({})".format(filter_name, df.shape[0]))
                continue

            for tool, *args in thresholds:
                try:
                    if np.sum(~df[tool].isnull()) == 0:
                        continue
                    generate_statistics(df, statistics, filter_name, tool, is_single_label=True)

                except KeyError:
                    continue

                stats_df = pd.DataFrame(statistics)
                stats_df.sort_values([metric], ascending=False).to_csv(outfile,
                                                                       sep="\t",
                                                                       index=False)
