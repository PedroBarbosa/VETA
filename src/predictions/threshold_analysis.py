from collections import OrderedDict, defaultdict
import pandas as pd
import numpy as np
from src.preprocessing.utils import ratio
import os
from src.plots.plots_threshold_analysis import plot_threshold
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

BETA_CONFIGS = [1, 2, 3, 5, 10]


def generate_threshold_analysis(dataset, filters, threshold_list, outdir, nofpoints=100):
    logging.info("----------------------------")
    logging.info("Starting threshold analysis based on Clinvar data")
    logging.info("----------------------------")
    outdir = os.path.join(outdir, "thresholds_analysis")
    os.mkdir(outdir)
    thresholds_to_return = {}
    for b in BETA_CONFIGS:
        thresholds_to_return[b] = {}
        for filtern, *args in filters:
            thresholds_to_return[b][filtern] = {}

    for filtern, filterf in filters:
        df_f = filterf(dataset).copy()
        if df_f.shape[0] < 60 or not all(i >= 30 for i in df_f['label'].value_counts().tolist()):
            logging.warning('WARN: {} Df has not minimum required size (and class balance) for threshold analysis'.format(filtern))
            continue

        final_thresholds = OrderedDict()
        for tool, direction, recommended_threshold, *args in threshold_list:
            if tool not in df_f.columns:
                continue
            df_ = df_f.loc[pd.notnull(df_f[tool]),].copy()
            if df_.shape[0] < 100:
                logging.warning("WARN: Not enough predicted variants (100) in the '{}' set by '{}' ({}). "
                                "Skipping thresholds analysis.".format(filtern, tool, df_.shape[0]))
                continue
            max_thr = df_[tool].max()
            min_thr = df_[tool].min()

            if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
                logging.warning("Something strange in max/min thresholds {} {} {}".format(tool, max_thr, min_thr))
                continue
            step = (max_thr - min_thr) / float(nofpoints)
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
                # df_ = df.copy()

                correct = np.sum(classification.eq(df_['label']))
                total = df_.shape[0]
                acc = ratio(correct, total)

                tp = np.sum(classification.eq(df_['label']) & classification)
                fp = np.sum(~df_['label'] & classification)
                fn = np.sum(df_['label'] & ~classification)
                tn = np.sum(classification.eq(df_['label']) & ~classification)
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

                # print(acc,sensitivity,precision,specificity)
                f1.append(ratio(2.0 * (precision * sensitivity), (sensitivity + precision)))

                r = ratio(np.sum(df_f['label']), df_f.shape[0])

            bf1, new_t = defaultdict(list), {}
            final_thresholds[tool] = []
            # print("New thresholds for {} tool".format(tool))
            for N in BETA_CONFIGS:
                for i, t in enumerate(threshold_range):
                    bf1[N].append(ratio((1.0 + N ** 2) * (precisions[i] * sensitivities[i]),
                                        (sensitivities[i] + (N ** 2) * precisions[i])))

                new_t[N] = threshold_range[bf1[N].index(max(bf1[N]))]
                # print("N={}:\t{}".format(N,new_t[N]))
                final_thresholds[tool].append(round(float(new_t[N]), 2))

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
                outdir
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
                outdir,
                simple=True

            )
        with open(os.path.join(outdir, "proposed_thresholds_{}.tsv".format(filtern)), 'w') as out:
            for tool, new_thresholds in final_thresholds.items():
                out.write("{}\t{}\n".format(tool, '\t'.join([str(t) for t in new_thresholds])))

        with open(os.path.join(outdir, "generated_proposed_thresholds_{}.tex".format(filtern)), 'w') as out:
            out.write("""\\begin{tabular}{ p{3cm} >{\\raggedleft\\arraybackslash}p{1.5cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm} >{\\raggedleft\\arraybackslash}p{1cm}}
\\hline
Tool       & Original & 1/1   & 1/2   & 1/3   & 1/5   & 1/10  \\\\
\\hline
""")
            for tool, new_thresholds in final_thresholds.items():
                original = ["{} {}".format(t[1], t[2]) for t in threshold_list if t[0] == tool][0]
                out.write(
                    "{}&\t{}&\t{}\n".format(tool, original, '\t&'.join([str(t) for t in new_thresholds])) + "\\\\\n")
            out.write("\\end{tabular}")

        for tool in final_thresholds:
            for beta, t in zip(BETA_CONFIGS, final_thresholds[tool]):
                thresholds_to_return[beta][filtern][tool] = t
    return thresholds_to_return
