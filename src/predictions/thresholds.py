import logging
import os
from collections import OrderedDict, defaultdict
from typing import final

import numpy as np
import pandas as pd

from plots.plots_threshold_analysis import plot_optimal_thresholds
from predictions.metrics import generate_statistics
from preprocessing.tables import generate_proposed_thresholds_latex
from preprocessing.osutils import ensure_folder_exists
from sklearn.metrics import f1_score


def perform_threshold_analysis(dataset: pd.DataFrame,
                               location_filters: list,
                               threshold_list: list,
                               tools_config: dict,
                               outdir: str,
                               f1_at_ref_thresh: dict,
                               n_of_points: int = 100,
                               beta_values: list = [1, 2, 3, 5]):
    """
    Do threshold analysis using Clinvar as ground
    truth dataset for threshold evaluation.

    Additionally, several precision/recall trade-offs
    are tested to find the most appropriate threshold
    per tool that maximimizes the number of true positives
    detected.

    :param pd.DataFrame dataset: Input dataset
    :param list location_filters: List of location filters
    :param list threshold_list: List of thresholds
    :param dict tools_config: Dict with config for each tool
    :param str outdir: Output directory
    :param dict f1_at_ref_thresh: Dict with F1 scores
    using the reference thresholds
    :param list beta_values: List of values to use as beta
        in the F_beta score calculation (precision/recall
        tradeoff). Default: `[1, 2, 3, 5]
    :param n_of_points: Default: `100`

    :return dict: New thresholds for all the tools at all locations
    """
    logging.info("---------------------------")
    logging.info("Starting threshold analysis")
    logging.info("---------------------------")
    outdir = os.path.join(outdir, "thresholds_analysis")
    ensure_folder_exists(outdir)
    
    available_for_threshold_analysis = ['coding', 'missense', 'splice_site', 'splice_region', 'intronic', 'deep_intronic']
  
    # Output dict initialization
    thresholds_to_return = {}
    for beta in beta_values:
        thresholds_to_return[beta] = {}
        for _loc in location_filters:
            if _loc not in available_for_threshold_analysis:
                continue
            thresholds_to_return[beta][_loc] = {}

    for _loc in location_filters:

        if _loc not in available_for_threshold_analysis:
            continue
        
        df_f = dataset[dataset.location == _loc].copy()
        print()
        logging.info("-------------------------")
        logging.info("Looking at {} variants (N={}, {} pos and {} neg)".format(_loc,
                                                                               df_f.shape[0],
                                                                               df_f[df_f.label].shape[0],
                                                                               df_f[~df_f.label].shape[0]))
        logging.info("-------------------------")
        # Do threshold analysis only if dataset has a minimum size
        if df_f.shape[0] < 60:
            logging.info('No minimum number of variants required (N=60)')
            continue
            
        if not all(i >= 30 for i in df_f.label.value_counts().tolist()):
            logging.warning('No minimum number of variants in the minority class required (N=30)')
            continue

        final_thresholds = OrderedDict()
        for tool, direction, reference_threshold, *args in threshold_list:

            # Threshold analysis is not performed for 
            # S-CAP (it uses multiple reference thresholds) and 
            # cVEP (predictions are discrete categories that VETA converts to fixed floats)
            if tool not in df_f.columns or tool.lower() in ['scap', 's-cap', 'cvep']:
                continue

            # Do threshold analysis for tools with sufficient prediction power
            df_per_tool = df_f.loc[pd.notnull(df_f[tool]),][[tool, 'label']].copy()

            # Transform TraP scores for proper threshold analysis
            if tool == "TraP":
                df_per_tool['TraP'] = df_per_tool['TraP'].apply(lambda x: x - 1 if x > 1 else x)

                # If all variants, it is not possible to distinguish coding/non-coding thresholds
                if _loc in ["coding", "missense"]:
                    reference_threshold = 0.416
                else:
                    reference_threshold = 0.289

            ratio_predicted = df_per_tool.shape[0] / df_f.shape[0]
            if df_per_tool.shape[0] < 60:
                logging.info('{} did not predict the minimum number of variants required (N=60)'.format(tool))
                continue
            
            if ratio_predicted < 0.5:
                logging.info("'{}' predicted less than 50% of variants (N={} predictions). "
                             "Threshold analysis will be performed, but be cautious in the "
                             "interpretation.".format(tool, df_per_tool.shape[0]))

            max_thr = df_per_tool[tool].max()
            min_thr = df_per_tool[tool].min()
 
            if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
                logging.info("Something strange in max/min thresholds {} {} "
                                "{}".format(tool, max_thr, min_thr))
                continue

            step = (max_thr - min_thr) / float(n_of_points)
            threshold_range = np.arange(min_thr, max_thr, step)

            statistics = defaultdict(list)
            # Compute metrics at each threshold
            for threshold in threshold_range:
                if direction == ">":
                    classification_f = lambda x: x == np.nan and np.nan or x >= threshold
                else:
                    classification_f = lambda x: x == np.nan and np.nan or x <= threshold

                df_per_tool[tool + '_prediction'] = df_per_tool[tool].map(classification_f)
                statistics = generate_statistics(df_per_tool, statistics, _loc, tool, f_beta=beta_values)

            new_t = {}
 
            final_thresholds[tool] = [(reference_threshold, f1_at_ref_thresh[_loc][tool])]

            accuracies = statistics['accuracy']
            sensitivities = statistics['sensitivity']
            specificities = statistics['specificity']
            f_betas = defaultdict(list)
            for _beta in beta_values:
                _res_at_each_t = statistics['Fbeta_' + str(_beta)]
                f_betas[_beta] = _res_at_each_t
                new_t[_beta] = threshold_range[_res_at_each_t.index(max(_res_at_each_t))]
                final_thresholds[tool].append((round(float(new_t[_beta]), 2), max(_res_at_each_t)))

            for v in [False, True]:
                plot_optimal_thresholds(tool,
                                        _loc,
                                        threshold_range,
                                        accuracies,
                                        sensitivities,
                                        specificities,
                                        f_betas,
                                        reference_threshold,
                                        new_t,
                                        outdir,
                                        simple=v)
        
        _outdir_new_thresh = os.path.join(outdir, "new_thresholds")
        _outdir_new_config = os.path.join(outdir, "new_configs")
        ensure_folder_exists(_outdir_new_config)
        ensure_folder_exists(_outdir_new_thresh)
        
        ###########################
        ### Write new thresholds ##
        ###########################
        generate_proposed_thresholds_latex(final_thresholds, _loc, _outdir_new_thresh)
        with open(os.path.join(_outdir_new_thresh, "proposed_thresholds_{}.tsv".format(_loc)), 'w') as out:
            out.write("Tool\tReference_threshold\tF1_score_at_ref_thresh\t{}\n".format("\t".join("1/{}".format(b) +
                                                                                                 "\t" + "Fbeta_{}".format(b)
                                                                                                 for b in beta_values)))
            for tool, new_thresholds in final_thresholds.items():
                out.write(
                    "{}\t{}\n".format(tool, '\t'.join(["{}\t{}".format(str(t[0]), str(t[1])) for t in new_thresholds])))

        for i, beta in enumerate(beta_values):

            config = open(os.path.join(_outdir_new_config, "new_config_{}_f{}.txt".format(_loc, beta)), 'w')

            for tool, new_thresh in final_thresholds.items():
                
                # delete original threshold
                _new_thresh = new_thresh[1:]
                thresholds_to_return[beta][_loc][tool] = _new_thresh
                
                t = _new_thresh[i][0]
                info = [v for v in threshold_list if v[0] == tool][0]
                vcf_field = ','.join([v for sublist in tools_config[tool] for v in sublist])
                outline = [tool, vcf_field, info[1], str(t), info[3]]
                config.write('\t'.join(outline) + '\n')
                
    return thresholds_to_return
