import os
from builtins import list

import matplotlib
import matplotlib.pyplot as plt
#matplotlib.use('agg')
import numpy as np
import matplotlib.patches as mpatches
from matplotlib import colors
from collections import OrderedDict, defaultdict
plt.switch_backend('agg')

from .utils import *

BETA_CONFIGS = [1, 2, 3, 5, 10]

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
    

def generate_threshold_analysis(dataset, filters, threshold_list, name, folder, nofpoints = 100):
    thresholds_to_return = {}
    for b in BETA_CONFIGS:
        thresholds_to_return[b] = {}
        for filtern, *args in filters:
            thresholds_to_return[b][filtern] = {}
    
    for filtern, filterf in filters:
        df_f = filterf(dataset).copy()
        if df_f.shape[0] < 60 or not all(i >= 30 for i in df_f['class'].value_counts().tolist()):
            print('{} Df has not minimum required size (and class balance) for threshold analysis'.format(filtern))
            continue
        set_style()
        final_thresholds = OrderedDict()
        for tool, direction, recommended_threshold, color, marker in threshold_list:
            if tool not in df_f.columns:
                continue
            df_ = df_f.loc[pd.notnull(df_f[tool]),].copy()
            if df_.shape[0] < 100:
                print("Not enough predicted variants (100) in the '{}' set by '{}' ({}). Skipping thresholds analysis.".format(filtern, tool, df_.shape[0]))
                continue
            max_thr = df_[tool].max()
            min_thr = df_[tool].min()

            if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
                print("Something strange in max/min thresholds {} {} {}".format(tool,max_thr,min_thr))
                continue
            step = (max_thr-min_thr)/float(nofpoints)
            threshold_range = np.arange(min_thr, max_thr, step)

            accuracies = []
            sensitivities = []
            specificities = []
            precisions = []
            f1 = []
            for threshold in threshold_range:
                if direction == ">":
                    classification_f = lambda x: x == np.nan and np.nan or x > threshold
                else:
                    classification_f = lambda x: x == np.nan and np.nan or x < threshold

                classification = df_[tool].map(classification_f)
                #df_ = df.copy()

                correct = np.sum(classification.eq(df_['class']))
                total = df_.shape[0]
                acc = ratio(correct, total)

                tp = np.sum(classification.eq(df_['class']) & classification)
                fp = np.sum(~df_['class'] & classification)
                fn = np.sum(df_['class'] & ~classification)
                tn = np.sum(classification.eq(df_['class']) & ~classification)
                ap = tp + fn
                ap_predicted = np.sum(classification)  # == tp +fp, sensitivity was being calculated with this value

                sensitivity = ratio(tp, ap)  # same as recall
                precision = ratio(tp, ap_predicted)

                an = tn + fp
                an_predicted = np.sum(~classification)  # == tn + fn, sensitivity was being calculated with this value
                specificity = ratio(tn, an)


                accuracies.append(acc)
                sensitivities.append(sensitivity)
                precisions.append(precision)
                specificities.append(specificity)

                #print(acc,sensitivity,precision,specificity)
                f1.append(ratio(2.0 * (precision * sensitivity), (sensitivity + precision)))

                r = ratio(np.sum(df_f['class']), df_f.shape[0])

            bf1, new_t = defaultdict(list),{}
            final_thresholds[tool] = []
            #print("New thresholds for {} tool".format(tool))
            for N in BETA_CONFIGS:
                for i, t in enumerate(threshold_range):
                    bf1[N].append(ratio((1.0 + N**2) * (precisions[i] * sensitivities[i]), (sensitivities[i] + (N**2)*precisions[i])))

                new_t[N] = threshold_range[bf1[N].index(max(bf1[N]))]
                #print("N={}:\t{}".format(N,new_t[N]))
                final_thresholds[tool].append(round(float(new_t[N]),2))

            plot_threshold(
                    tool,
                    filtern,
                    threshold_range,
                    accuracies,
                    sensitivities,
                    specificities,
                    f1,
                    bf1,
                    recommended_threshold,
                    new_t,
                    folder
            )
            plot_threshold(
                    tool,
                    filtern,
                    threshold_range,
                    accuracies,
                    sensitivities,
                    specificities,
                    f1,
                    bf1,
                    recommended_threshold,
                    new_t,
                    folder,
                    simple = True
                
            )
        with open(os.path.join(folder, 'datasets', "proposed_thresholds_{}.tsv".format(filtern)),'w') as out:
            for tool, new_thresholds in final_thresholds.items():
                out.write("{}\t{}\n".format(tool,'\t'.join([str(t) for t in new_thresholds])))
                
                
        with open(os.path.join(folder, "generated_proposed_thresholds_{}.tex".format(filtern)),'w') as out:
            out.write("""\\begin{tabular}{ p{3cm} >{\\raggedleft\\arraybackslash}p{1.5cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm}}
\\hline
Tool       & Original & 1/1   & 1/2   & 1/3   & 1/5   & 1/10  \\\\
\\hline
""")
            for tool, new_thresholds in final_thresholds.items():
                
                original = [ "{} {}".format(t[1],t[2]) for t in threshold_list if t[0] == tool][0]
                out.write("{}&\t{}&\t{}\n".format(tool, original, '\t&'.join([str(t) for t in new_thresholds])) + "\\\\\n")
            out.write("\\end{tabular}")                
        
        for tool in final_thresholds:
            for beta, t in zip(BETA_CONFIGS, final_thresholds[tool]):
                thresholds_to_return[beta][filtern][tool] = t
    return thresholds_to_return
    


