import seaborn as sns
sns.set(style="white")
cmap = sns.diverging_palette(220, 10, as_cmap=True)
from plots.plots_performance_comparison import *
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


def classify_all_variants(df, thresholds):
    for tool, direction, threshold, color, marker in thresholds:
        if direction == ">":
            classification = lambda x: pd.isna(x) and np.nan or x > threshold
        else:
            classification = lambda x: pd.isna(x) and np.nan or x < threshold

        prediction = df[tool].apply(classification)
        df[tool + '_prediction'] = prediction
    return df


def apply_new_thresholds(dataset, filtername, thresholds, new_thresholds):
    if new_thresholds:
        ntlines = []
        for t in thresholds:
            ntline = [x for x in t]

            if t[0] in new_thresholds[filtername]:
                ntline[2] = new_thresholds[filtername][t[0]]
            ntlines.append(ntline)

        df = classify_all_variants(dataset, ntlines)
        suffix = "_proposed"
    else:
        suffix = "_original"
        df = classify_all_variants(dataset, thresholds)
    return df, suffix


def generate_performance_comparison(dataset, filtes_var_type, filters, thresholds, name, folder, new_thresholds=None):
    if new_thresholds:
        logging.info("#################################")
        logging.info("Tools performance analysis for {} dataset with the new thresholds started".format(name))
        logging.info("#################################")
    else:
        logging.info("#################################")
        logging.info("Tools performance analysis for {} dataset started".format(name))
        logging.info("#################################")

    for vartype, vartypefunction in filtes_var_type:

        ensure_folder_exists(os.path.join(folder, "figures", vartype))
        outdir = os.path.join(folder, "figures", vartype)
        df_v = vartypefunction(dataset).copy()
        logging.info("Looking at {} ({} variants)".format(vartype, df_v.shape[0]))
        if df_v.shape[0] == 0:
            logging.warning("WARN: There are no {} in the variant set. Skipping this analysis.".format(vartype))
            continue

        for filtername, filterfunction in filters:
            statistics = defaultdict(list)
            df = filterfunction(df_v).copy()
            if df.shape[0] < 10:
                logging.warning("WARN: {} has not a minimum number of {} {} variants (10) to evaluate tools performance."
                                " ({})".format(name, vartype, filtername, df.shape[0]))
                continue
            df, suffix = apply_new_thresholds(df, filtername, thresholds, new_thresholds)

            for tool, *args in thresholds:
                if (vartype == "snps" or vartype == "all_types") and (filtername == "all" or filtername == "splicesite"):
                    plot_density_by_class(df[[tool, "class"]], thresholds, os.path.join(outdir, 'tools_analysis_' +
                                                                                        name + "_" + filtername + "_" +
                                                                                        tool + suffix))
                s_df = df[~df[tool + '_prediction'].isnull()]
                if np.sum(~df[tool + '_prediction'].isnull()) == 0:
                    continue
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
                statistics['accuracy'].append(accuracy)
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

            stats_df = pd.DataFrame(statistics)
            plot_tools(stats_df, os.path.join(outdir, 'tools_analysis_' + name + "_" + filtername + suffix))
            plot_tools_paper(stats_df, os.path.join(outdir, 'tools_analysis_paper_' + name + "_" + filtername + suffix))
            plot_metrics(stats_df,
                         os.path.join(outdir, 'tools_performance_metrics_' + name + "_" + filtername + suffix))
            stats_df.drop(['filter'], axis=1).to_csv(os.path.join(outdir, "statistics_{}.csv").format(filtername),
                                                     sep="\t", index=False)
            if "weighted_accuracy" in list(stats_df):
                stats_df[["tool", "weighted_accuracy", "accuracy", "weighted_f1", "f1"]].sort_values(
                    ["weighted_accuracy"], ascending=False).to_csv(
                    os.path.join(outdir, "tools_ranking_{}.csv").format(filtername), sep="\t", index=False)
            else:
                stats_df[["tool", "weighted_f1", "f1", "accuracy"]].sort_values(["weighted_f1"],
                                                                                ascending=False).to_csv(
                    os.path.join(outdir, "tools_ranking_{}.csv").format(filtername), sep="\t", index=False)

    logging.info("Done!")
