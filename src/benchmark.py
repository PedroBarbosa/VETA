import fnmatch
from collections import defaultdict
from typing import List
import pandas as pd

from base import Base
from plots.plots_benchmark_mode import *
from plots.plots_machine_learning import plot_feature_correlation
from predictions import introns
from predictions import metrics
from predictions.apply import apply_tool_predictions
from predictions.metapredictions import Metapredictions
from predictions.thresholds import perform_threshold_analysis
from preprocessing.clinvar import *
from preprocessing.tables import generate_clinvar_table, generate_consequence_table
from preprocessing.osutils import ensure_folder_exists


class BenchmarkTools(Base):
    """
    Benchmark mode to compare tools
    based on a reference dataset
    """

    def __init__(self, dataset: str,
                 out_dir: str,
                 scopes_to_predict: List = None,
                 types_of_variant: List = None,
                 metric: str = "weighted_accuracy",
                 location: str = "HGVSc",
                 aggregate_classes: bool = False,
                 select_conseqs: str = "gene_body",
                 do_intronic_analysis: bool = False,
                 split_splice_sites: bool = False,
                 clinvar_stars: str = "1s_l",
                 phenotype_ids: list = None,
                 do_threshold_analysis: bool = False,
                 do_bootstrapping: bool = False,
                 do_machine_learning: bool = False,
                 allele_frequency_col: str = "gnomADg_AF",
                 skip_heatmap: bool = False,
                 tools_config: str = None):

        """
        ----
        Base args described in Base class
        ----
        """

        dataset, is_clinvar = self.check_dataset_arg(dataset)
        super().__init__(vcf=dataset,
                         out_dir=out_dir,
                         scopes_to_predict=scopes_to_predict,
                         types_of_variant=types_of_variant,
                         metric=metric,
                         location=location,
                         aggregate_classes=aggregate_classes,
                         select_conseqs= select_conseqs,
                         do_intronic_analysis=do_intronic_analysis,
                         split_splice_sites=split_splice_sites,
                         is_clinvar=is_clinvar,
                         allele_frequency_col=allele_frequency_col,
                         skip_heatmap=skip_heatmap,
                         tools_config=tools_config)

        self.clinvar_stars = clinvar_stars
        self.do_threshold_analysis = do_threshold_analysis
        self.do_bootstrapping = do_bootstrapping
        self.do_machine_learning = do_machine_learning

        # Remove tools where all values are missing
        tools = [t[0] for t in self.thresholds]
        self.df.dropna(axis=1, how="all", inplace=True)
        self.null_tools = list(set(tools) - set(list(self.df)))
        if len(self.null_tools) > 0:
            logging.info("Tools with no predictions were found. They will be "
                         "discarded from analysis: {}".format(self.null_tools))

        if self.is_clinvar:
            self.out_dir = os.path.join(self.out_dir, self.clinvar_stars)

            self.df = filter_by_condition(self.df, phenotype_ids)
            
            df_clinvar_sure = filter_clinvar_sure(self.df)
            self.clinvar_levels_dict = {
                '0s': df_clinvar_sure,
                '1s': filter_clinvar_1_star(df_clinvar_sure),
                '2s': filter_clinvar_2_stars(df_clinvar_sure),
                '3s': filter_clinvar_3_stars(df_clinvar_sure),
                '4s': filter_clinvar_4_stars(df_clinvar_sure),
                '0s_l': self.df,
                '1s_l': filter_clinvar_1_star(self.df),
                '2s_l': filter_clinvar_2_stars(self.df),
                '3s_l': filter_clinvar_3_stars(self.df),
                '4s_l': filter_clinvar_4_stars(self.df)}

            generate_clinvar_table(self.clinvar_levels_dict, 
                                   self.variant_types, 
                                   self.out_dir,
                                   self.clinvar_stars)

            self.df = self.clinvar_levels_dict[self.clinvar_stars]
            assert self.df.shape[0] > 0, "No remaining variants after filtering by {} " \
                                         "level. Try a more relaxed value for the '-c'" \
                                         " arg.".format(self.clinvar_stars)
            logging.info('Number of variants after filtering by Clinvar stars ({}): {}'.format(self.clinvar_stars, self.df.shape[0]))

        generate_consequence_table(self.df, 
                                   self.out_dir)

        self.top_tools, f1_at_ref_threshold = self.do_performance_comparison()
        if self.do_intronic_analysis:
            thresholds = [tool for tool in self.thresholds if tool[3] != 'Protein']
         
            introns.IntronicAnalysis(self.df,
                                     self.split_splice_sites,
                                     thresholds,
                                     self.metric,
                                     self.aggregate_classes,
                                     self.out_dir,
                                     self.allele_frequency_col)

        if self.do_threshold_analysis:
            
            new_thresholds = perform_threshold_analysis(self.df,
                                                        self.location_filters,
                                                        self.thresholds,
                                                        self.tools_config,
                                                        self.out_dir,
                                                        f1_at_ref_threshold,
                                                        do_bootstrapping=self.do_bootstrapping)

        if self.do_machine_learning:
            self.do_ml_analysis()

    def do_performance_comparison(self):
        """
        Generates performance stats and plots
        for the processed df

        :return dict: df with predictions
        """
        logging.info("----------------------------------")
        logging.info("Tools performance analysis started")
        logging.info("----------------------------------")

        top_tools, f1_at_ref_threshold = {}, {}
        if self.is_clinvar:
            logging.info("Evaluations are done with Clinvar variants that belong "
                         "to the filtering strategy employed ({}). Playing with '--clinvar_stars' "
                         "argument allows using more permissive or restrictive variant sets.".format(self.clinvar_stars))
            
            
        for var_type, _var_type_func in self.variant_types:
            outdir = os.path.join(self.out_dir, "tools_benchmark", var_type)
            ensure_folder_exists(outdir)

            _df_v = _var_type_func(self.df).copy()
            print()
            logging.info("----------")
            logging.info("Looking at {} ({} variants)".format(var_type, _df_v.shape[0]))
            logging.info("----------")
            print()
            
            if _df_v.shape[0] == 0:
                logging.warning("No variants. Skipping!".format(var_type))
                continue
            
            for _location in self.location_filters:
                
                if _location == "all":
                    df = _df_v.copy()
                else:
                    df = _df_v[_df_v.location == _location].copy()
                    
                if df.empty:
                    continue
                
                statistics = defaultdict(list)
                        
                n_pos = df.label.sum()
                n_neg = (~ df.label).sum()
                df['count_class'] = df.groupby('outcome')['outcome'].transform('size')
                
                logging.info("----------")
                logging.info("{} variants (n={})".format(_location, df.shape[0]))
                logging.info("Pathogenic: {}".format(n_pos))
                logging.info("Benign: {}".format(n_neg))
                logging.info("----------")
                print()
                
                if df.shape[0] < 10:
                    logging.info("Less than 10 '{}' '{}' variants found ({})"
                                ". Skipping!".format(var_type, _location, df.shape[0]))
                    continue
                
                df = self.apply_thresholds(df, _location)

                if not self.skip_heatmap:
                    if df.shape[0] > 10:
                        logging.info("Generating performance heatmaps for '{}' '{}' variants. "
                                    "If dataset is large, this step will take long time.".format(var_type, _location))
                        
                        fname = os.path.join(outdir, "performance_at_fixed_thresh", "heatmap_" + _location + '.pdf')
                        plot_heatmap(df, fname=fname)
                        logging.info("Done")
                    else:
                        logging.info("Less than 10 variants located in '{}' of type '{}'. " \
                                    "Skipping heatmap generation.".format(_location, var_type))

                ensure_folder_exists(os.path.join(outdir, "results_tsv"))
                out_preds = os.path.join(outdir, "results_tsv", "preds_{}_{}_class.tsv").format(var_type, _location)
                out_preds_scores = os.path.join(outdir, "results_tsv", "preds_{}_{}_scores.tsv").format(var_type, _location)
                out_c = ['chr', 'pos', 'ref', 'alt', 'SYMBOL'] + [x for x in list(df) if '_prediction' in x] + ['label']
                out_c_scores = ['chr', 'pos', 'ref', 'alt', 'SYMBOL'] + [x for x in list(df) if '_prediction' not in x and x not in ['count_class', 'blank', 'label']]
  
                df[out_c].to_csv(out_preds, index=False, sep="\t")
                df[out_c_scores].to_csv(out_preds_scores, index=False, sep="\t")
                roc_metrics_per_tool, pr_metrics_per_tool = [], []

                for tool, direction, _, _ in self.thresholds:

                    try:

                        statistics = metrics.generate_statistics(df, statistics, _location, tool)
                        scored = df[~df[tool + '_prediction'].isnull()]
                        if scored.empty:
                            logging.info("{} did not score any {} variant.".format(tool, _location))
                            continue
                       
                        if 0 < scored.shape[0] < df.shape[0] * 0.05:
                            logging.info("Warn: {} scored some '{}' variants (scored={}), but they account for less than "
                                        "5% of the dataset size. Plots will not show these results.".format(tool, _location, scored.shape[0]))
                            continue

                        # Distributions of each class are only plotted for SNPs and
                        # all types for splice_site and coding locations
                        if (_location in ['all', 'splice_site', 'splice_region', 'intronic', 'missense', 'coding'] and var_type in ['snps', 'all_types']):
                            plot_density_by_class(df[[tool, "label"]],
                                                thresholds=self.thresholds,
                                                fname=os.path.join(outdir, 'class_distribution', tool + "_" + _location))

                        ####################
                        ### ROC analysis ###
                        ####################
                        # Required to have predictions on variants of both classes
                        n_pos_pred = scored[scored.label].shape[0]
                        n_neg_pred = scored[scored.label == False].shape[0]

                        if scored.empty or 0 < scored.shape[0] < df.shape[0] * 0.33:
                            logging.info("Warn: {} scored some '{}' variants (scored={}), but they account for less than "
                                        "33% of the dataset size. ROC analysis will be skipped for this tool.".format(tool, _location, scored.shape[0]))
                        
                        elif min(n_pos_pred, n_neg_pred) < 15:
                            logging.info("Warn: No minimum number of predictions on each class (15) found. Skipping ROC analysis for {} tool in {} variant set".format(tool, _location))
                            
                        else:
                    
                            na_frac = 1 - (scored.shape[0] / df.shape[0])

                            try:
                                roc_curve, pr_curve, roc_auc, ap_score = metrics.do_roc_analysis(scored[[tool, 'label']],
                                                                                                tool,
                                                                                                higher_is_better=direction == ">")
            
                                na_frac = 1 - (scored.shape[0] / df.shape[0])
                                roc_metrics_per_tool.append([tool, na_frac, roc_curve[0], roc_curve[1],
                                                            roc_curve[2], roc_curve[3], roc_auc])
                                pr_metrics_per_tool.append([tool, na_frac, pr_curve[0], pr_curve[1],
                                                            pr_curve[2], pr_curve[3], ap_score])
                            except TypeError:

                                roc_metrics_per_tool.append([tool, na_frac, None, None, None, None, None])
                                pr_metrics_per_tool.append([tool, na_frac, None, None, None, None, None])
                                pass
                        
                    except KeyError:
                        pass
                
                stats_all_df = pd.DataFrame(statistics)
    
                if roc_metrics_per_tool:
                    roc_m = pd.DataFrame([[x[0], x[-1]] for x in roc_metrics_per_tool], columns=['tool', 'auROC'])
                    pr_m = pd.DataFrame([[x[0], x[-1]] for x in pr_metrics_per_tool], columns=['tool', 'average_precision'])
                    stats_all_df = stats_all_df.merge(roc_m, on='tool', how='left').merge(pr_m, on='tool', how='left')

                # draw plots
                out_stats = os.path.join(outdir, "results_tsv", "statistics_{}_{}.tsv").format(var_type, _location)
                out_ranks = os.path.join(outdir, "results_tsv", "tools_ranking_{}_{}.tsv").format(var_type, _location)

                af_plot = os.path.join(outdir, "allele_frequency", "AF_{}_{}.pdf".format(var_type, _location))
                unscored_plot = os.path.join(outdir, "performance_at_fixed_thresh", "unscored_fraction_{}_{}.pdf".format(var_type, _location))
                barplot_plot = os.path.join(outdir, "performance_at_fixed_thresh", "barplot_{}_{}.pdf".format(var_type, _location))
                metrics_plot = os.path.join(outdir, "performance_at_fixed_thresh", "scatter_{}_{}.pdf".format(var_type, _location))

                plot_allele_frequency(df, af_plot, self.allele_frequency_col)
                plot_unscored(stats_all_df, unscored_plot)
                plot_tools_barplot(stats_all_df, barplot_plot, self.metric)
                plot_metrics(stats_all_df, metrics_plot, self.metric)

                for i, _roc_data in enumerate([roc_metrics_per_tool, pr_metrics_per_tool]):

                    if i == 0:
                        is_roc = True 
                        bname = "ROC_"
                    else:
                        is_roc = False
                        bname = "ROC_pr_"
                    
                    if _roc_data:
                        
                        _ = plot_curve(_roc_data, 
                                os.path.join(outdir, "roc_analysis", "{}{}.pdf".format(bname, _location)),
                                (n_pos, n_neg),
                                is_roc=is_roc)

                stats_all_df.drop(['filter'], axis=1).to_csv(out_stats, sep="\t", index=False)

                ranked_stats = stats_all_df[["tool", "weighted_accuracy", "accuracy",
                                            "weighted_F1", "F1", "coverage",
                                            "specificity", "sensitivity",
                                            "precision", "mcc", "norm_mcc", "weighted_norm_mcc"]].sort_values([self.metric], ascending=False)

                top_tools[_location] = ranked_stats.head(20) if ranked_stats.shape[0] > 20 else ranked_stats
                f1_at_ref_threshold[_location] = dict(zip(ranked_stats.tool, ranked_stats.F1))
                ranked_stats.to_csv(out_ranks, sep="\t", index=False)

        logging.info("Done!")
        return top_tools, f1_at_ref_threshold

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

    def do_ml_analysis(self):
        """
        Machine learning analysis of tools scores.
        """
        logging.info("---------------------------------")
        logging.info("Starting meta-prediction analysis")
        logging.info("---------------------------------")
        out_dir = os.path.join(self.out_dir, "machine_learning")
        ensure_folder_exists(out_dir)

        for _loc in self.location_filters:
            print()
            logging.info("----------")
            logging.info("Analyzing '{}' variants.".format(_loc))
            logging.info("----------")
            if _loc == "all":
                df_f = self.df.copy()
            else:
                df_f = self.df[self.df.location == _loc].copy()

            if df_f.shape[0] < 100:
                logging.info("Less than 100 variants in '{}' dataframe (N={}). "
                             "Machine learning analysis will be skipped.".format(_loc, df_f.shape[0]))
                continue
            
            n_pos_pred = df_f[df_f.label].shape[0]
            n_neg_pred = df_f[df_f.label == False].shape[0]
            if min(n_pos_pred, n_neg_pred) < 30:
                logging.info("Less than 30 variants of the minority class (only {}) in '{}' "
                             "variants. Classifiers will not be trained.".format(min(n_pos_pred, n_neg_pred), _loc))
                continue
            
            ml_data = Metapredictions(df_f, _loc, self.thresholds,
                                      out_dir,
                                      top_tools=self.top_tools,
                                      rank_metric=self.metric)

            plot_feature_correlation(ml_data.df, _loc, out_dir)

            if len(ml_data.features) < 3:
                logging.info("Less than 3 features (aka tool scores) available"
                             " ({}) after removing predictors with many missing "
                             "values in '{}' variants. Classifiers will not be "
                             "trained for such small array size.".format(len(ml_data.features),
                                                                         _loc))
                continue

            logging.info("Applying classifiers to data.")
            trained_clfs = ml_data.do_classification()

            logging.info("Generating feature ranking.")
            ml_data.generate_feature_ranking(trained_clfs)
            logging.info("Generating decision tree plot")
            ml_data.generate_tree_plot(trained_clfs)
