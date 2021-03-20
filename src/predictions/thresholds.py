import logging
import os
from collections import OrderedDict, defaultdict

import numpy as np
import pandas as pd

from src.plots.plots_threshold_analysis import plot_optimal_thresholds
from src.predictions.metrics import generate_statistics
from src.preprocessing.latex import generate_proposed_thresholds_latex
from src.preprocessing.osutils import ensure_folder_exists


def perform_threshold_analysis(dataset: pd.DataFrame,
                               location_filters: list,
                               threshold_list: list,
                               outdir: str,
                               n_of_points: int = 100,
                               beta_values: list = [1, 2, 3, 5, 10]):
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
    :param str outdir: Output directory
    :param list beta_values: List of values to use as beta
        in the F_beta score calculation (precision/recall
        tradeoff). Default: `[1, 2, 3, 5, 10]
    :param n_of_points: Default: `100`

    :return dict: New thresholds for all the tools at all locations
    """
    logging.info("---------------------------")
    logging.info("Starting threshold analysis")
    logging.info("---------------------------")
    outdir = os.path.join(outdir, "thresholds_analysis")
    ensure_folder_exists(outdir)

    # Output dict initialization
    thresholds_to_return = {}
    for beta in beta_values:
        thresholds_to_return[beta] = {}
        for location, *args in location_filters:
            thresholds_to_return[beta][location] = {}

    for _loc, filter_func in location_filters:

        df_f = filter_func(dataset).copy()

        # Do threshold analysis only if dataset has a minimum size
        if df_f.shape[0] < 60 or not all(i >= 30 for i in df_f['label'].value_counts().tolist()):
            logging.warning('WARN: Dataframe of \'{}\' variants does not have the minimum '
                            'required size (60) and/or class balance (> 30 variants of '
                            'minority class) for threshold analysis'.format(_loc))
            continue

        final_thresholds = OrderedDict()
        for tool, direction, reference_threshold, *args in threshold_list:

            # Threshold analysis is not performed for S-CAP (it uses multiple reference thresholds)
            if tool not in df_f.columns or tool == "SCAP" or tool == "S-CAP":
                continue

            # Do threshold analysis for tools with sufficient prediction power
            df_per_tool = df_f.loc[pd.notnull(df_f[tool]), ][[tool, 'label']].copy()

            # Transform TraP scores for proper threshold analysis
            if tool == "TraP":
                df_per_tool['TraP'] = df_per_tool['TraP'].apply(lambda x: x - 1 if x > 1 else x)

                # If all variants, it is not possible to distinguish coding/non-coding thresholds
                if _loc == "all":
                    continue
                if _loc == "coding":
                    reference_threshold = 0.416
                else:
                    reference_threshold = 0.289

            ratio_predicted = df_per_tool.shape[0] / df_f.shape[0]
            if ratio_predicted < 0.5 and df_per_tool.shape[0] < 60:
                logging.warning("WARN: '{}' predicted less than 50% of '{}' variants "
                                "({} predictions). Skipping thresholds analysis because "
                                "less than 60 variants were predicted.".format(tool, _loc, df_per_tool.shape[0]))
                continue

            elif ratio_predicted < 0.5:
                logging.warning("WARN: '{}' predicted less than 50% of '{}' variants "
                                "({} predictions). Thresholds analysis will be performed, "
                                "but be cautious in the interpretation.".format(tool, _loc, df_per_tool.shape[0]))

            max_thr = df_per_tool[tool].max()
            min_thr = df_per_tool[tool].min()

            if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
                logging.warning("Something strange in max/min thresholds {} {} "
                                "{}".format(tool, max_thr, min_thr))
                continue

            step = (max_thr - min_thr) / float(n_of_points)
            threshold_range = np.arange(min_thr, max_thr, step)

            statistics = defaultdict(list)
            # Compute metrics at each threshold
            for threshold in threshold_range:
                if direction == ">":
                    classification_f = lambda x: x == np.nan and np.nan or x > threshold
                else:
                    classification_f = lambda x: x == np.nan and np.nan or x < threshold

                df_per_tool[tool + '_prediction'] = df_per_tool[tool].map(classification_f)
                statistics = generate_statistics(df_per_tool, statistics, _loc, tool, f_beta=beta_values)

            new_t = {}
            final_thresholds[tool] = [reference_threshold]
            accuracies = statistics['accuracy']
            sensitivities = statistics['sensitivity']
            specificities = statistics['specificity']
            f_betas = defaultdict(list)
            for _beta in beta_values:
                _res_at_each_t = statistics['Fbeta_' + str(_beta)]
                f_betas[_beta] = _res_at_each_t
                new_t[_beta] = threshold_range[_res_at_each_t.index(max(_res_at_each_t))]
                final_thresholds[tool].append(round(float(new_t[_beta]), 2))

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

        generate_proposed_thresholds_latex(final_thresholds, _loc, outdir)
        with open(os.path.join(outdir, "proposed_thresholds_{}.tsv".format(_loc)), 'w') as out:
            out.write("Tool\tReference_threshold\t{}\n".format("\t".join("1/" + str(b) for b in beta_values)))
            for tool, new_thresholds in final_thresholds.items():
                out.write("{}\t{}\n".format(tool, '\t'.join([str(t) for t in new_thresholds])))

        for tool in final_thresholds:
            # slicing = delete original threshold
            for beta, t in zip(beta_values, final_thresholds[tool][1:]):
                thresholds_to_return[beta][_loc][tool] = t
    return thresholds_to_return
