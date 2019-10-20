import seaborn as sns

sns.set(style="white")
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from plots.plots_performance_comparison import *
from plots.plots_intronic_analysis import *
import sys
import logging
from sklearn.metrics import auc
from scipy import integrate

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(message)s')


def classify_all_variants(df, thresholds):
    for tool, direction, threshold, *args in thresholds:
        if direction == ">":
            classification = lambda x: pd.isna(x) and np.nan or x > threshold
        else:
            classification = lambda x: pd.isna(x) and np.nan or x < threshold

        try:
            prediction = df[tool].apply(classification)
            df[tool + '_prediction'] = prediction
        except KeyError:
            continue
    return df


def generate_statistics(df, statistics, filtername, tool):
    s_df = df[~df[tool + '_prediction'].isnull()]
    statistics['filter'].append(filtername)
    statistics['tool'].append(tool)

    tp = np.sum(s_df['class'].eq(s_df[tool + '_prediction']) & s_df['class'])
    tn = np.sum(s_df['class'].eq(s_df[tool + '_prediction']) & ~s_df['class'])

    fp = np.sum(s_df['class'].ne(s_df[tool + '_prediction']) & ~s_df['class'])
    fn = np.sum(s_df['class'].ne(s_df[tool + '_prediction']) & s_df['class'])

    mp = np.sum(df[tool + '_prediction'].isnull() & df['class'])
    mn = np.sum(df[tool + '_prediction'].isnull() & ~df['class'])
    nan = np.sum(df[tool + '_prediction'].isnull())

    total = df.shape[0]
    correct = tp + tn
    precision = ratio(tp, tp + fp)
    recall = ratio(tp, tp + fn)
    accuracy = ratio(correct, (total - nan))
    coverage = ratio(tp + tn + fp + fn, total)

    statistics['total'].append(total)
    statistics['correct'].append(correct)
    statistics['nan'].append(nan)
    statistics['fraction_nan'].append(ratio(nan, total))
    statistics['coverage'].append(coverage)
    statistics['accuracy'].append(round(accuracy, 2))
    statistics['precision'].append(precision)
    statistics['specificity'].append(ratio(tn, tn + fp))
    statistics['sensitivity'].append(recall)
    statistics['tp'].append(tp)
    statistics['fp'].append(fp)
    statistics['tn'].append(tn)
    statistics['fn'].append(fn)
    statistics['mp'].append(mp)
    statistics['mn'].append(mn)

    statistics['scored_p'].append(np.sum(s_df['class'] == True))
    statistics['scored_n'].append(np.sum(s_df['class'] == False))

    statistics['total_p'].append(np.sum(df['class'] == True))
    statistics['total_n'].append(np.sum(df['class'] == False))

    try:
        statistics['f1'].append(round(2 * (precision * recall) / (precision + recall), 2))
        statistics['weighted_f1'].append(
            round((2 * (precision * recall) / (precision + recall)) * coverage, 2))
    except ZeroDivisionError:
        statistics['f1'].append(0)
        statistics['weighted_f1'].append(0)
    # if np.sum(df['class']) == np.sum(~df['class']):
    statistics['weighted_accuracy'].append(round(accuracy * coverage, 2))
    return statistics


def perform_intron_analysis(df, filter_intronic_bins, threshold_list, metric_to_evaluate, dataset_name, out_dir):
    print("\n")
    logging.info("---------------------------------------------")
    logging.info("Starting intronic variants analysis for {} dataset".format(dataset_name))
    logging.info("---------------------------------------------")
    os.mkdir(os.path.join(out_dir, "intron_analysis"))
    df_i = df[~df['intron_bin'].isnull()]
    booleanDictionary = {True: 'Pathogenic', False: 'Benign'}
    df_i["outcome"] = df_i["class"].map(booleanDictionary)
    plot_general_bin_info(df_i, filter_intronic_bins, os.path.join(os.path.join(out_dir, "intron_analysis"),
                                                                   "intronic".format(dataset_name)))

    roc_per_bin = defaultdict(list)
    na = {}
    for bin, filterbin in filter_intronic_bins:
        df_i_ = filterbin(df_i).copy()
        df_i_['count_class'] = df_i_.groupby('outcome')['outcome'].transform('size')

        if df_i_.shape[0] > 20:
            print("\n")
            logging.info("Looking at {} bin ({} variants)".format(bin, df_i_.shape[0]))
        else:
            logging.info("Not enough variants for ROC analysis at {} bin ({})".format(bin, df_i_.shape[0]))
            continue

        logging.info("Pathogenic: {}".format(df_i_["class"].sum()))
        logging.info("Benign: {}".format((~ df_i_["class"]).sum()))

        plot_allele_frequency(df_i_, os.path.join(os.path.join(out_dir, "intron_analysis"),
                                                  "AF_{}".format(bin)))
        list_df_metrics_per_tool = []
        statistics = defaultdict(list)

        for tool, direction, threshold, *args in threshold_list:
            try:
                df_ = df_i_.loc[pd.notnull(df_i_[tool]),].copy()
            except KeyError:
                na[tool] = 1
                continue

            statistics = generate_statistics(df_i_, statistics, bin, tool)
            stats_df = pd.DataFrame(statistics)
            na[tool] = stats_df.loc[stats_df['tool'] == tool, 'fraction_nan'].iloc[0]
            df_tool = df_[~df_[tool + "_prediction"].isnull()]
            if df_tool.shape[0] >= 0:

                max_thr = df_tool[tool].max()  # + (df_tool[tool].max()) - df_tool[tool].min()) * 0.001
                min_thr = df_tool[tool].min()  # - (df_tool[tool].max()) - df_tool[tool].min()) * 0.001

                if pd.isnull(max_thr) and pd.isnull(min_thr):
                    logging.info("{} has no predicted variants for {} bin".format(tool, bin))
                    continue
                elif max_thr == min_thr:
                    logging.info("{} has the the same min/max score: {}".format(tool, max_thr))
                    continue

                step = (max_thr - min_thr) / float(100)
                threshold_range = np.arange(min_thr, max_thr, step)

                tool_metrics = []
                for threshold in threshold_range:
                    if direction == ">":
                        classification_f = lambda x: x == np.nan and np.nan or x >= threshold
                    else:
                        classification_f = lambda x: x == np.nan and np.nan or x <= threshold

                    classification = df_[tool].map(classification_f)
                    tp = np.sum(classification.eq(df_['class']) & classification)
                    fp = np.sum(~df_['class'] & classification)
                    fn = np.sum(df_['class'] & ~classification)
                    tn = np.sum(classification.eq(df_['class']) & ~classification)
                    ap = tp + fn
                    ap_predicted = np.sum(classification)  # == tp +fp, sensitivity was being calculated with this value

                    sensitivity = ratio(tp, ap)
                    precision = ratio(tp, ap_predicted)

                    an = tn + fp
                    an_predicted = np.sum(
                        ~classification)  # == tn + fn, sensitivity was being calculated with this value
                    specificity = ratio(tn, an)
                    f1 = ratio(2.0 * (precision * sensitivity), (sensitivity + precision))

                    tool_metrics.append([tool + "(n=" + str(df_tool.shape[0]) + ",", threshold, precision, sensitivity,
                                         1 - specificity])

                prc = [val_at_thresh[2] for val_at_thresh in tool_metrics]
                fpr = [val_at_thresh[4] for val_at_thresh in tool_metrics]
                tpr = [val_at_thresh[3] for val_at_thresh in tool_metrics]

                roc_auc = [auc(sorted(fpr), sorted(tpr))] * len(tool_metrics)
                pr_auc = [auc(sorted(tpr), sorted(prc))] * len(tool_metrics)
                # pr_auc = [integrate.simps(prc, dx=step)] * len(tool_metrics)
                # trapz gives the same are, but with negative values. Need to reverse X
                # roc_auc = [np.trapz(tpr, fpr, dx=step)] * len(tool_metrics)
                # pr_auc = [np.trapz(prc, tpr)] * len(tool_metrics)

                f1_score = [stats_df.loc[stats_df['tool'] == tool, 'f1'].iloc[0]] * len(tool_metrics)
                weighted_f1_score = [stats_df.loc[stats_df['tool'] == tool, 'weighted_f1'].iloc[0]] * len(tool_metrics)
                nan = [na[tool]] * len(tool_metrics)
                tool_metrics = [x + [y] + [z] + [f] + [wf] + [na] for x, y, z, f, wf, na in zip(tool_metrics, roc_auc,
                                                                                                pr_auc,
                                                                                                f1_score,
                                                                                                weighted_f1_score,
                                                                                                nan)]
                list_df_metrics_per_tool.append(pd.DataFrame(tool_metrics,
                                                             columns=["tool", "threshold", "precision", "recall", "FPR",
                                                                      "ROC-auc", "PR-auc",
                                                                      "F1", "weighted_F1", "fraction_nan"]))

                if bin != "all_intronic" and bin is not None:
                    roc_per_bin[tool].append([bin, roc_auc[0], pr_auc[0], f1_score[0], weighted_f1_score[0], na[tool]])
            else:
                logging.info("\t{} didn't score at least 20 variants in the {} bin".format(tool, bin))

        stats_df.drop(['filter'], axis=1).to_csv(os.path.join(os.path.join(out_dir, "intron_analysis"),
                                                              "statistics_{}.csv".format(bin)), sep="\t", index=False)

        plot_unscored(stats_df,
                      os.path.join(os.path.join(out_dir, "intron_analysis"), 'unscored_fraction_{}'.format(bin)))
        plot_metrics(stats_df, os.path.join(os.path.join(out_dir, "intron_analysis"), 'performance_{}'.format(bin)),
                     metric_to_evaluate)
        final_df_metrics = pd.concat(list_df_metrics_per_tool)

        try:
            plot_ROCs(final_df_metrics, os.path.join(os.path.join(out_dir, "intron_analysis"), "{}_ROC".format(bin)),
                      df_i_.loc[df_i_['outcome'] == "Pathogenic", 'count_class'].iloc[0])
        except IndexError:
            logging.info("No positive class instances, skipping ROC")
    df_roc_bins = pd.DataFrame([[k] + i for k, v in roc_per_bin.items() for i in v], columns=["tool", "bin", "auROC",
                                                                                              "prROC",
                                                                                              "F1",
                                                                                              "weighted_F1",
                                                                                              "fraction_nan"])
    plot_auROC_by_bin(df_roc_bins, os.path.join(os.path.join(out_dir, "intron_analysis", "per_bin_evol")))


def apply_thresholds(dataset, filtername, thresholds, new_thresholds=None):

    if new_thresholds:
        ntlines = []
        for t in thresholds:
            ntline = [x for x in t]
            if t[0] in new_thresholds[filtername]:
                ntline[2] = new_thresholds[filtername][t[0]]
            ntlines.append(ntline)
        df = classify_all_variants(dataset, ntlines)
        suffix = "_proposed"
        return df, suffix

    else:
        df = classify_all_variants(dataset, thresholds)
        return df


def generate_performance_comparison(dataset, filtes_var_type, filters, thresholds, metric_to_evaluate, name, folder):
    logging.info("---------------------------------")
    logging.info("Tools performance analysis for {} dataset started".format(name))
    logging.info("---------------------------------")

    for vartype, vartypefunction in filtes_var_type:

        ensure_folder_exists(os.path.join(folder, "tools_benchmark", vartype))
        outdir = os.path.join(folder, "tools_benchmark", vartype)
        df_v = vartypefunction(dataset).copy()
        logging.info("Looking at {} ({} variants)".format(vartype, df_v.shape[0]))
        if df_v.shape[0] == 0:
            logging.warning("WARN: There are no {} in the variant set. Skipping this analysis.".format(vartype))
            continue

        for filtername, filterfunction in filters:
            statistics = defaultdict(list)
            df = filterfunction(df_v).copy()
            if df.shape[0] < 10:
                logging.info("{} has not a minimum number of {} {} variants (10) to evaluate tools performance. ({})".
                             format(name, vartype, filtername, df.shape[0]))
                continue

            df = apply_thresholds(df, filtername, thresholds)
            for tool, *args in thresholds:
                try:
                    if (vartype == "snps" or vartype == "all_types") and \
                            (filtername == 'all' or filtername == "splicesite"):
                        plot_density_by_class(df[[tool, "class"]],
                                              thresholds,
                                              os.path.join(outdir, 'class_dist_' + tool + "_" + filtername))

                    scored = np.sum(~df[tool + '_prediction'].isnull())
                    if 0 < scored < df.shape[0] / 10:
                        logging.info("{} scored some {} variants, but they account for less than "
                                     "10% of the dataset size. Although stats file will include these results,"
                                     "plots will not show those.".format(tool, filtername))

                    statistics = generate_statistics(df, statistics, filtername, tool)

                except KeyError:
                    pass

            stats_all_df = pd.DataFrame(statistics)
            stats_df = stats_all_df[stats_all_df.fraction_nan < 0.90]

            plot_unscored(stats_all_df,
                              os.path.join(outdir, 'unscored_fraction_' + filtername))
            plot_tools(stats_df, os.path.join(outdir, 'tools_analysis_' + filtername ), metric_to_evaluate)
            plot_tools_paper(stats_df, os.path.join(outdir, 'tools_analysis_paper_' + filtername))
            plot_metrics(stats_df,
                         os.path.join(outdir, 'tools_performance_metrics_' + filtername),
                         metric_to_evaluate)

            stats_all_df.drop(['filter'], axis=1).to_csv(os.path.join(outdir, "statistics_{}.csv").format(filtername),
                                                         sep="\t", index=False)

            stats_all_df[["tool", "weighted_accuracy", "accuracy", "weighted_f1", "f1", "coverage"]].sort_values(
                [metric_to_evaluate], ascending=False).to_csv(
                os.path.join(outdir, "tools_ranking_{}.csv").format(filtername), sep="\t", index=False)

    logging.info("Done!")


def generate_performance_comparison_with_new_thresholds(dataset, filters_var_type, filters, thresholds,
                                                        metric_to_evaluate, name, folder, new_thresholds):
    logging.info("---------------------------------")
    logging.info("Tools performance analysis for {} dataset with the new thresholds started".format(name))
    logging.info("---------------------------------")

    for vartype, vartypefunction in filters_var_type:

        ensure_folder_exists(os.path.join(folder, "thresholds_analysis", vartype))
        outdir = os.path.join(folder, "thresholds_analysis", vartype)
        df_v = vartypefunction(dataset).copy()
        logging.info("Looking at {} ({}) variants with new thresholds.".format(vartype, df_v.shape[0]))
        if df_v.shape[0] == 0:
            logging.warning("WARN: There are no {} in the variant set. Skipping this analysis.".format(vartype))
            continue

        filter_with_new_thresholds = [filter for filter in filters if filter[0] in ["all", "coding", "splicesite"]]
        for filtername, filterfunction in filter_with_new_thresholds:
            statistics = defaultdict(list)
            df = filterfunction(df_v).copy()
            if df.shape[0] < 10:
                logging.info("{} has not a minimum number of {} {} variants (10) to evaluate tools performance. ({})".
                             format(name, vartype, filtername, df.shape[0]))
                continue

            df, suffix = apply_thresholds(df, filtername, thresholds, new_thresholds)
            for tool, *args in thresholds:
                try:
                    statistics = generate_statistics(df, statistics, filtername, tool)
                except KeyError:
                    pass

            stats_all_df = pd.DataFrame(statistics)
            stats_df = stats_all_df[stats_all_df.fraction_nan < 0.90]

            plot_tools(stats_df, os.path.join(outdir, 'tools_analysis_' + filtername + suffix), metric_to_evaluate)
            plot_metrics(stats_df,
                         os.path.join(outdir, 'tools_performance_metrics_' + filtername + suffix),
                         metric_to_evaluate)

            stats_all_df.drop(['filter'], axis=1).to_csv(os.path.join(outdir, "statistics_{}_{}.csv").
                                                         format(filtername, suffix),
                                                         sep="\t", index=False)

            stats_all_df[["tool", "weighted_accuracy", "accuracy", "weighted_f1", "f1", "coverage"]].sort_values(
                [metric_to_evaluate], ascending=False).to_csv(
                os.path.join(outdir, "tools_ranking_{}_{}.csv").format(filtername, suffix),
                sep="\t", index=False)

    logging.info("Done!")