import os

import matplotlib.patches as mpatches
import numpy as np

from .plots_utils import *

plt.switch_backend('agg')


def plot_optimal_thresholds(tool: str, location: str, threshold_range: np.array,
                            accuracies: list, sensitivities: list, specificities: list,
                            beta_values: dict,
                            reference_threshold: float,
                            new_thresholds: dict,
                            outdir: str,
                            simple=False):
    """
    :param str tool: Tool to plot new thresholds
    :param str location: Location filter for the variants analysed
    :param np.array threshold_range: Observed range of scores for the tool
    :param list accuracies: Accuracy at each threshold value
    :param list sensitivities: Sensitivity/recall at each threshold value
    :param list specificities: Specificity at each threshold value
    :param dict beta_values: Fbeta scores at each threshold value for each
        beta values used in the analysis
    :param float reference_threshold: Reference threshold for the tool
    :param dict new_thresholds: Optimized threshold at each beta value
    :param str outdir: Output directory
    :param boolean simple: Draw simpler plot. Default: `False`. All new
        recommended thresholds at each beta value will be plotted.
    """

    plt.figure(figsize=(8, 6))
    plt.plot(threshold_range, accuracies, label="Accuracy", linewidth=1, color='#B22400')
    plt.plot(threshold_range, sensitivities, label="Sensitivity", linewidth=1, color='#006BB2')
    plt.plot(threshold_range, specificities, label="Specificity", linewidth=1, color='#429C40')
    plt.axvline(x=reference_threshold, color='b', linestyle=':', linewidth=2, label="Recommended threshold")

    patches = []
    if not simple:
        col = sns.color_palette("Greys_r", len(beta_values))
        col_i = 0
        for k in sorted(beta_values.keys()):
            plt.plot(threshold_range, beta_values[k], color=col[col_i], linestyle="--", linewidth=1)
            plt.axvline(x=new_thresholds[k], color=col[col_i], linestyle="--", linewidth=1)
            patches.append(mpatches.Patch(color=col[col_i], label='+1 TP = +{} FP'.format(k), linestyle="--"))
            col_i += 1

    plt.subplots_adjust(right=0.7)
    leg1 = plt.legend(bbox_to_anchor=(1, 0.7), frameon=True)

    if not simple:
        plt.legend(handles=patches, title="New thresholds", bbox_to_anchor=(1, 1))
    plt.gca().add_artist(leg1)
    plt.title("{}".format(tool))
    plt.xlabel("Threshold values")
    plt.ylabel("Metrics")

    f_name = os.path.join(outdir, "threshold_{}_{}{}.pdf".format(tool.replace(" ", "_"),
                                                                 location,
                                                                 simple and "_simple" or ""))
    plt.savefig(f_name, bbox_inches='tight')
    plt.gcf().clear()
    plt.clf()
    plt.cla()
    plt.close()
