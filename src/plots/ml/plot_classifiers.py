import os.path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')
sns.set(style="white")
from plots.plots_utils import *
import matplotlib.patches as mpatches


def plot_classifiers(out_metrics, filtern, folder, extra=None):
    fig, axes = plt.subplots(ncols=2, sharey=True, figsize=(15, 15))
    classifier_names = [k for k in out_metrics.keys() if np.average([p['accuracy'] for p in out_metrics[k]]) > 0]
    ind = np.arange(len(classifier_names))

    correct_p = [np.average([p['tp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    incorrect_p = [np.average([p['fn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    correct_n = [np.average([p['tn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    incorrect_n = [np.average([p['fp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]

    missing_p = [np.average([p['mp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    missing_n = [np.average([p['mn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]

    w = 0.8

    axes[0].barh(ind, correct_p, align='center', color='darkblue', zorder=10, height=w, linewidth=0, alpha=0.7)
    axes[0].barh(ind, incorrect_p, left=correct_p, align='center', color='darkred', zorder=10, height=w, linewidth=0, alpha=0.7)
    axes[0].barh(ind, missing_p, left=[x + y for (x, y) in zip(incorrect_p, correct_p)], align='center',
                 color='lightgrey', zorder=10, height=w, linewidth=0, alpha=0.9)

    axes[1].barh(ind, correct_n, align='center', color='darkblue', zorder=10, height=w, linewidth=0, alpha=0.7)
    axes[1].barh(ind, incorrect_n, left=correct_n, align='center', color='darkred', zorder=10, height=w, linewidth=0, alpha=0.7)
    axes[1].barh(ind, missing_n, left=[x + y for (x, y) in zip(incorrect_n, correct_n)], align='center',
                 color='lightgrey', zorder=10, height=w, linewidth=0, alpha=0.9)

    axes[0].invert_xaxis()
    axes[0].set(yticks=ind,
                yticklabels=["{} ({:.1f}%)".format(c, np.average([p['accuracy'] for p in out_metrics[c]]) * 100) for c
                             in classifier_names])

    axes[0].set_xlabel('Pathogenic Variants')
    axes[1].set_xlabel('Benign Variants')
    for ax in axes.flat:
        ax.grid(False)
        ax.margins(0.00)
        ax.tick_params(axis='both', which='major', labelsize=8)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)

    colors = ["darkblue", "darkred", "lightgrey"]
    patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
    axes[1].legend(patches, ["Correct", "Incorrect", "Unpredictable"], bbox_to_anchor=(1.5, 1),
                   borderaxespad=0, loc='upper right', prop=dict(size=8))
    set_size(fig, len(classifier_names))
    fig.tight_layout()
    if extra:
        extra = "_" + extra
    else:
        extra = ""
    plt.savefig(os.path.join(folder,  "classifiers_comparison_{}{}.pdf".format(filtern, extra)))
    plt.close()
