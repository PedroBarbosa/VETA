import logging
import pandas as pd
from collections import defaultdict
from typing import List
from preprocessing.osutils import ensure_folder_exists
from plots.plots_intronic_analysis import *
from plots.plots_benchmark_mode import plot_allele_frequency, plot_unscored, plot_metrics, plot_curve
from predictions import metrics
from predictions.apply import apply_tool_predictions


def do_intron_analysis(df: pd.DataFrame, 
                       thresholds: List,
                       metric: str,
                       aggregate_classes: str,
                       out_dir: str,
                       af_column: str,
                       min_variants: int = 20):
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
    assert "intron_bin" in df.columns, "Intronic bins not in the data. " \
                                       "Probably a first run of Clinvar " \
                                       "was performed without --intronic-bins " \
                                       "args. To fix, just remove the 'tsv' in " \
                                       "the input directory and try again."

    out_dir = os.path.join(out_dir, "intron_analysis")
    out_tsv = os.path.join(out_dir, 'results_tsv')
    out_af = os.path.join(out_dir, "allele_frequency")
    out_fixed_thresh = os.path.join(out_dir, "performance_at_fixed_thresh")
    out_roc = os.path.join(out_dir, "roc_analysis")
    out_all_bin_agg = os.path.join(out_dir, "all_bin_together")
    [ensure_folder_exists(path) for path in [out_dir, out_tsv, out_af, out_fixed_thresh, out_roc, out_all_bin_agg]]
 
    # intronic variants
    df_i = df[~df['intron_bin'].isnull()].copy(deep=True)
    df_i = apply_tool_predictions(df_i, thresholds)

    roc_per_bin = defaultdict(list)
    na = {}
    loc_in_protein_coding = ['splice_site', 'splice_region', 'intronic', 'deep_intronic']
    pc_is_done = False
    bin_to_exclude = ['1-2', '3-10'] if aggregate_classes else ['1-10']
    locations = df_i.location.unique().tolist()
    locations.append('all')

    for _class in locations:
    
        if _class == "all":
            _df_loc = df_i.copy(deep=True)
        elif _class in loc_in_protein_coding:
            if pc_is_done is False:
                _class = "protein_coding_regions"
                _df_loc = df_i[df_i.location.isin(loc_in_protein_coding)].copy(deep=True)
                pc_is_done = True
            else:
                continue
        else:
            _df_loc = df_i[df_i.location == _class].copy(deep=True)
        
        if _class != "protein_coding_regions":
            continue
        
        print("\n\n\n")
        logging.info("--------------")
        logging.info("Looking at {} variants".format(_class))
        logging.info("--------------")    
        print()
        if _df_loc.empty:
            logging.info('No {} variants found.'.format(_class))
            continue
        
        _df_loc['count_class'] = _df_loc.groupby('outcome')['outcome'].transform('size')
        
        try:
            plot_general_bin_info(_df_loc, 
                                    out_dir, 
                                    _class, 
                                    aggregate_classes, 
                                    af_column)
        except ValueError:
            logging.info("Problem plotting info about intronic bins. Skipping.")

        for _bin, _filter_func in filter_intronic_bins:
            if "all_except" in _bin or _bin == "all_intronic":
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
            n_neg = (~ _df_i_bin["label"]).sum()
            logging.info("Pathogenic: {}".format(n_pos))
            logging.info("Benign: {}".format(n_neg))
                
            if _df_i_bin.shape[0] < min_variants:
                logging.info("Not enough {} variants for ROC analysis " 
                             "at {} bin (N={}, required={})".format(_class, _bin, _df_i_bin.shape[0], min_variants))
                skip_roc = True
            
            elif min(n_pos, n_neg) < 10:
                logging.info("Number of variants in the minority class is lower (N={}) "
                            "than minimum required (N=10) for ROC analysis.".format(min(n_pos, n_neg)))
                skip_roc = True
                    
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
                if skip_roc is False:
                    n_pos_pred = df_tool[df_tool[tool + "_prediction"]].shape[0]
                    n_neg_pred = df_tool[df_tool[tool + "_prediction"] == False].shape[0]

                    if float(na[tool]) > 0.5:
                        logging.info("ROC analysis will be skipped for {} at {} bin. "
                                    "More than 50% of missing data ({})".format(tool, _bin, round(float(na[tool]), 2)))

                    elif min(n_pos_pred, n_neg_pred) < 10:
                        logging.info("ROC analysis will be skipped for {} at {} bin. "
                                    "No minimum number of predictions on each class (10) found (N={})".format(tool, _bin, min(n_pos_pred, n_neg_pred)))
                        
                    else:

                        try:
                            roc_curve, pr_curve, roc_auc, ap_score = metrics.do_roc_analysis(df_tool[[tool, 'label']],
                                                                                            tool)
                            roc_metrics_per_tool.append([tool, float(na[tool]), roc_curve[0], roc_curve[1],
                                                        roc_curve[2], roc_curve[3], roc_auc])
                            pr_metrics_per_tool.append([tool, float(na[tool]), pr_curve[0], pr_curve[1],
                                                        pr_curve[2], pr_curve[3], ap_score])
                        # S-cap
                        except TypeError:
                            pass


                f1_score = stats_df.loc[stats_df['tool'] == tool, 'F1'].iloc[0]
                weighted_f1_score = stats_df.loc[stats_df['tool'] == tool, 'weighted_F1'].iloc[0]
                if not _bin in ["all_intronic", "all_except_1-2", "all_except_1-10", None]:
                    roc_per_bin[tool].append([_class,
                                              _bin,
                                              f1_score,
                                              weighted_f1_score,
                                              na[tool]])

                stats_out = os.path.join(out_tsv, "statistics_{}_{}.tsv".format(_class, _bin))
                stats_df.drop(['filter'], axis=1).to_csv(stats_out, sep="\t", index=False)
                af_plot = os.path.join(out_af, "AF_{}_{}.pdf".format(_class, _bin))
                unscored_plot = os.path.join(out_fixed_thresh, "unscored_fraction_{}_{}.pdf".format(_class, _bin))
                metrics_plot = os.path.join(out_fixed_thresh, "tools_metrics_{}_{}.pdf".format(_class, _bin))

                plot_allele_frequency(_df_i_bin, af_plot, af_column)
                plot_unscored(stats_df, unscored_plot)
                plot_metrics(stats_df, metrics_plot, metric)

                if roc_metrics_per_tool:
                    plot_curve(roc_metrics_per_tool,
                            os.path.join(out_roc, "{}_{}_ROC.pdf".format(_class, _bin)),
                            (n_pos, n_neg))

                if pr_metrics_per_tool:
                    plot_curve(pr_metrics_per_tool,
                            os.path.join(out_roc, "{}_{}_ROC_pr.pdf".format(_class, _bin)),
                            (n_pos, n_neg),
                            is_roc=False)

        df_roc_bins = pd.DataFrame([[k] + i for k, v in roc_per_bin.items() for i in v],
                                   columns=["tool", "variant_class", "bin", "F1", "weighted_F1", "fraction_nan"])
        
        df_roc_bins.groupby('variant_class').apply(plot_metrics_by_bin, 
                                                   os.path.join(out_all_bin_agg, "per_bin"), 
                                                   aggregate_classes)
