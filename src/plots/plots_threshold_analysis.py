import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
plt.switch_backend('agg')
from .plots_utils import *

def plot_threshold(
            tool,
            filtern,
            thresholds,
            accuracies,
            sensitivities,
            specificities,
            f1,
            bf1,
            recommended_threshold,
            new_threshold,
            outdir,
            simple=False
        ):

    plt.figure(figsize=(8,6))
    plt.plot(thresholds, accuracies, label="Accuracy", linewidth=1, color='#B22400')
    plt.plot(thresholds, sensitivities, label="Sensitivity", linewidth=1,  color='#006BB2')
    plt.plot(thresholds, specificities, label="Specificity", linewidth=1,  color='#429C40')
    #plt.plot(thresholds, f1, label="F1", color='darkviolet')
    plt.axvline(x=recommended_threshold, color='b',linestyle=':', label="Recommended threshold")
    
    if not simple:
        patches=[]
        col = sns.color_palette("Greys_r", len(bf1))
        #col[-1] = colors.to_rgba('antiquewhite', alpha=None)
        col_i=0
        for k in sorted(bf1.keys()):
            plt.plot(thresholds, bf1[k], color=col[col_i], linestyle="--", linewidth=1)
            plt.axvline(x=new_threshold[k], color=col[col_i],linestyle="--", linewidth=1)
            patches.append(mpatches.Patch(color=col[col_i],label='+1 TP = +{} FP'.format(k), linestyle="--"))
            col_i+=1

    plt.subplots_adjust(right=0.7)
    leg1=plt.legend(loc="best",frameon=True)
    if not simple:
        plt.legend(handles=patches,title="New thresholds",bbox_to_anchor=(1, 0.5))
    plt.gca().add_artist(leg1)
    plt.title("Threshold analysis for {}".format(tool))
    plt.xlabel("Threshold values")
    plt.ylabel("Metrics")
    fname = os.path.join(outdir, 'figures', "threshold_analysis_{}_{}{}.pdf".format(tool.replace(" ", "_"), filtern, simple and "_simple" or ""))
    plt.savefig(fname, bbox_inches='tight')
    plt.gcf().clear()
    plt.clf()
    plt.cla()
    plt.close()
    


    


