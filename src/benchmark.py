from typing import List
from src.predictions.filters import filters_location
from src.base import Base
from src.predictions.apply import apply_tool_predictions
from src.preprocessing.osutils import ensure_folder_exists
from src.plots.plots_benchmark_mode import *
import pandas as pd
from src.predictions import introns
from src.predictions import metrics
from collections import defaultdict
import sys
import logging
import fnmatch

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')
from src.preprocessing.latex import generate_clinvar_table
from src.preprocessing.clinvar import *
TOOLS_CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "map_tools2vcf_annotation.txt")


class BenchmarkTools(Base):
    """
    Benchmark mode to compare tools
    based on a reference dataset
    """

    def __init__(self, dataset: str,
                 out_dir: str,
                 scope_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 genome: str = "hg19",
                 is_intronic: bool = False,
                 clinvar_stars: str = "3s_l",
                 do_threshold_analysis: bool = False,
                 do_machine_learning: bool = False,
                 skip_heatmap: bool = False,
                 tools_config: str = TOOLS_CONFIG):

        """
        ----
        Base args described in Base class
        ----

        :param  str clinvar_stars:

        :param bool do_threshold_analysis:

        :param bool do_machine_learning:

        :param bool skip_heatmap:
        """
        dataset, is_clinvar = self.check_dataset_arg(dataset)
        super().__init__(vcf=dataset,
                         out_dir=out_dir,
                         scope_to_predict=scope_to_predict,
                         types_of_variant=types_of_variant,
                         metric=metric,
                         location=location,
                         genome=genome,
                         is_intronic=is_intronic,
                         is_clinvar=is_clinvar,
                         tools_config=tools_config)

        self.clinvar_stars = clinvar_stars
        self.do_threshold_analysis = do_threshold_analysis
        self.do_machine_learning = do_machine_learning

        # Remove tools with where all values are missing
        tools = [t[0] for t in self.thresholds]
        self.df.dropna(axis=1, how="all", inplace=True)
        self.null_tools = list(set(tools) - set(list(self.df)))
        if len(self.null_tools) > 0:
            logging.info("Tools with no predictions were found. They will be "
                         "discarded from analysis: {}".format(self.null_tools))

        if self.is_clinvar:
            self.out_dir = os.path.join(self.out_dir, self.clinvar_stars)
            df_clinvar_sure = filter_clinvar_sure(self.df)

            self.clinvar_levels_dict = {
                'clinvar': df_clinvar_sure,
                '1s': filter_clinvar_1_star(df_clinvar_sure),
                '2s': filter_clinvar_2_stars(df_clinvar_sure),
                '3s': filter_clinvar_3_stars(df_clinvar_sure),
                '4s': filter_clinvar_4_stars(df_clinvar_sure),
                'clinvar_l': self.df,
                '1s_l': filter_clinvar_1_star(self.df),
                '2s_l': filter_clinvar_2_stars(self.df),
                '3s_l': filter_clinvar_3_stars(self.df),
                '4s_l': filter_clinvar_4_stars(self.df)}

            generate_clinvar_table(self.clinvar_levels_dict, self.variant_types, self.out_dir)

            self.df = self.clinvar_levels_dict[self.clinvar_stars]
            assert self.df.shape[0] > 0, "No remaining variants after filtering by {} " \
                                         "level. Try a more relaxed value for the '-c'" \
                                         " arg.".format(self.clinvar_stars)

        self.do_performance_comparison()

        if self.is_intronic:
            introns.do_intron_analysis(self.df,
                                       thresholds=self.thresholds,
                                       metric=self.metric,
                                       out_dir=self.out_dir)

    def do_performance_comparison(self):
        """
        Generates performance stats and plots
        for the processed df

        :param pd.DataFrame df: Input df
        :param dict kwargs: Arbitrary number
            of args with all the attributes from
            of the ToolBenchmark object
        :return pd.DataFrame: df with predictions
        """
        logging.info("----------------------------------")
        logging.info("Tools performance analysis started")
        logging.info("----------------------------------")

        for var_type, _var_type_func in self.variant_types:
            outdir = os.path.join(self.out_dir, "tools_benchmark", var_type)
            ensure_folder_exists(outdir)

            _df_v = _var_type_func(self.df).copy()
            logging.info("Looking at {} ({} variants)".format(var_type, _df_v.shape[0]))
            if _df_v.shape[0] == 0:
                logging.warning("No variants. Skipping!".format(var_type))
                continue

            for _location, _filter_function in self.location_filters:
                statistics = defaultdict(list)
                df = _filter_function(_df_v).copy()
                df['count_class'] = df.groupby('outcome')['outcome'].transform('size')

                if df.shape[0] < 10:
                    logging.info("Less than 10 {} {} variants found ({})"
                                 ". Skipping!".format(var_type, _location, df.shape[0]))
                    continue

                df = self.apply_thresholds(df, _location)

                for tool, *args in self.thresholds:
                    try:
                        # distributions of each class
                        # are only plotted for SNPs and
                        # all types for splicesite and
                        # coding locations
                        if (filters_location == 'all' or
                            filters_location == "splicesite" or
                            filters_location == "coding") and (var_type == "snps" or
                                                               var_type == "all_types"):

                            plot_density_by_class(df[[tool, "label"]],
                                                  thresholds=self.thresholds,
                                                  fname=os.path.join(outdir,
                                                                     'class_dist_' + tool
                                                                     + "_" + _location))

                        scored = np.sum(~df[tool + '_prediction'].isnull())
                        if 0 < scored < df.shape[0] / 10:
                            logging.info("{} scored some {} variants, but they account for less than "
                                         "10% of the dataset size. Although stats file will include these results,"
                                         "plots will not show those.".format(tool, _location))

                        statistics = metrics.generate_statistics(df, statistics, _location, tool)

                    except KeyError:
                        pass

                stats_all_df = pd.DataFrame(statistics)
                stats_df = stats_all_df[stats_all_df.fraction_nan < 0.90]

                # draw plots
                out_stats = os.path.join(outdir, "statistics_{}_{}.tsv").format(var_type, _location)
                out_ranks = os.path.join(outdir, "tools_ranking_{}_{}.csv").format(var_type, _location)
                af_plot = os.path.join(outdir, "AF_{}_{}".format(var_type, _location))
                unscored_plot = os.path.join(outdir, "unscored_fraction_{}_{}".format(var_type, _location))
                barplot_plot = os.path.join(outdir, "tools_analysis_{}_{}".format(var_type, _location))
                metrics_plot = os.path.join(outdir, "tools_metrics_{}_{}".format(var_type, _location))

                plot_allele_frequency(df, af_plot)
                plot_unscored(stats_all_df, unscored_plot)
                plot_tools_barplot(stats_df, barplot_plot, self.metric)
                plot_metrics(stats_df, metrics_plot, self.metric)

                stats_all_df.drop(['filter'], axis=1).to_csv(out_stats, sep="\t", index=False)
                stats_all_df[["tool", "weighted_accuracy", "accuracy",
                              "weighted_F1", "F1", "coverage"]].sort_values([self.metric],
                                                                            ascending=False).to_csv(out_ranks,
                                                                                                    sep="\t",
                                                                                                    index=False)
        logging.info("Done!")

    def apply_thresholds(self, df: pd.DataFrame,
                         filter_location: str,
                         new_thresholds: pd.DataFrame = None):
        """
        Apply tools predictions based on a
        list of thresholds

        :param pd.DataFrame df: Input df
        :param str filter_location: Location filter
        :param pd.DataFrame new_thresholds: Df of new
            thresholds to evaluate predictions.
            Only applicable after a threshold
            analysis was performed by setting
            `--do_threshold_analysis` in VETA.

        :return pd.DataFrame: Df with tools predictions
            based on a list of reference thresholds
        """

        if new_thresholds is not None:
            ntlines = []

            for _t in new_thresholds:
                ntline = [x for x in _t]

                if _t[0] in new_thresholds[filter_location]:
                    ntline[2] = new_thresholds[filter_location][_t[0]]
                ntlines.append(ntline)

            df = apply_tool_predictions(df, ntlines)
            suffix = "_proposed"
            return df, suffix

        else:
            return apply_tool_predictions(df, self.thresholds)

    def check_dataset_arg(self, dataset: str):
        """
        Checks the input dataset argument
        :param str dataset:
        :return str: Dataset info
        :return bool: Whether dataset
            refers to clinvar (if a file)
        """
        # if reference datasets are in directory
        if os.path.isdir(dataset):
            fname_benign = [filename for filename in os.listdir(dataset) if
                            any(fnmatch.fnmatch(filename, pattern) for pattern in ["*benign*[vcf.bgz]",
                                                                                   "*neutral*[vcf.bgz]"])]
            fname_patho = [filename for filename in os.listdir(dataset) if
                           any(fnmatch.fnmatch(filename, pattern) for pattern in ["*deleterious*[vcf.bgz]",
                                                                                  "*pathogenic*[vcf.bgz]"])]
            if fname_benign and fname_patho:
                benign = os.path.join(dataset, fname_benign[0])
                pathogenic = os.path.join(dataset, fname_patho[0])
                return (benign, pathogenic), False
            else:
                raise ValueError('Files with required patterns absent in {}'.format(dataset))

        # if a file, veta expects to be from Clinvar
        elif os.path.isfile(dataset):
            return dataset, True

        else:
            raise ValueError("Set a valid value for the input dataset "
                             "(directory or file)")
