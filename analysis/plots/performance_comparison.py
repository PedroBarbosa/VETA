import numpy as np
import seaborn as sns

sns.set(style="white")
from .utils import *
import os
from collections import defaultdict
from osutils import ensure_folder_exists
from matplotlib import lines

# matplotlib.use('agg')
plt.switch_backend('agg')
cmap = sns.diverging_palette(220, 10, as_cmap=True)


def plot_density_by_class(data, thresholds, fname):
    fraction_nan = round(np.sum(np.sum(data.iloc[:, 0].isnull())) / data.shape[0] * 100, 2)
    df_ = data[~data.iloc[:, 0].isnull()]
    if df_.shape[0] > 0:
        pathogenic = df_.loc[df_['class'] == True]
        benign = df_.loc[df_['class'] == False]
        sns.distplot(benign.iloc[:, 0], hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3}, label="Benign")
        sns.distplot(pathogenic.iloc[:, 0], hist=False, kde=True,
                     kde_kws={'shade': True, 'linewidth': 3}, label="Pathogenic")
        plt.axvline([x[2] for x in thresholds if x[0] == list(df_)[0]], color='r', linestyle="--")

        # thresh_line = lines.Line2D([], [], color='r', linestyle='--',
        #                             markersize=10, markeredgewidth=1.5, label='Reference threshold')
        # plt.legend(handles=[thresh_line])

        plt.xlabel(list(df_)[0] + " ({}% missing data)".format(fraction_nan))
        ends = (max(df_.iloc[:, 0]) - min(df_.iloc[:, 0])) * 0.05
        plt.xlim(min(df_.iloc[:, 0]) - ends, max(df_.iloc[:, 0]) + ends)
        plt.savefig(fname + ".pdf")
        plt.close()


def plot_metrics(data, fname):
    my_range = range(1, len(data.index) + 1)
    data.sort_values('weighted_accuracy', ascending=True, inplace=True)
    fig, ax = plt.subplots(figsize=(10, 10))
    #  plt.hlines(y=my_range, xmin=data['specificity'], xmax=data['sensitivity'], color='grey', alpha=0.75)
    plt.scatter(data['specificity'], my_range, color='skyblue', alpha=1, marker='s', edgecolors='black', linewidths=0.5,
                label='Specificity')
    plt.scatter(data['sensitivity'], my_range, color='yellow', alpha=0.75, marker='^', edgecolors='black',
                linewidths=0.5, label='Sensitivity')
    plt.scatter(data['coverage'], my_range, color='brown', alpha=0.75, marker='o', edgecolors='black', linewidths=0.5,
                label='Fraction_predictable')

    i = 0
    for idx, rows in data.iterrows():
        val_range = list(rows[['specificity', 'sensitivity', 'coverage']])
        width = [min(val_range), max(val_range)]
        height = 1
        ax.add_patch(
            plt.Rectangle(xy=(width[0] - 0.01, my_range[i] - 0.4), width=width[1] - width[0] + 0.02, height=height,
                          linewidth=0.75, color='gray', fill=False))
        i += 1

    ax.grid(axis=0, linestyle='dashed')
    plt.legend()
    plt.yticks(my_range, data['tool'] + " (" + data['weighted_accuracy'].astype(str) + ")")
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools_by_type(data, fname):
    df = data.iloc[::-1]
    set_style()
    ind = np.arange(df.shape[0])
    df.sort_values(by=['correct'], inplace=True)
    w = 0.8

    fig, ax = plt.subplots()
    plt.xlim(df['total'][0])
    ax.barh(ind, df['correct'], align='center', color='green', zorder=10, alpha=0.8, height=w, edgecolor="black",
            linewidth=0.75)
    ax.set(yticks=ind, yticklabels=df['tool'])
    ax.grid(True)
    ax.margins(0.00)
    ax.invert_xaxis()
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xlabel('# Number of correctly predicted variants')
    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools_paper(data, fname):
    df = data.iloc[::-1]
    set_style()

    w = 0.5
    fig, axes = plt.subplots(ncols=2, sharey=True)
    df.sort_values('tp', inplace=True, ascending=True)
    ind = np.arange(df.shape[0])

    axes[0].barh(ind, df['tp'], align='center', color='green', zorder=10, height=w, edgecolor='black', linewidth=0.75)
    axes[1].barh(ind, df['tn'], align='center', color='green', zorder=10, height=w, edgecolor='black', linewidth=0.75)

    axes[0].set_xlim(0, df["tp"][0] + df["fn"][0])
    axes[1].set_xlim(0, df["tn"][0] + df["fp"][0])
    axes[0].set(yticks=ind, yticklabels=df['tool'])
    axes[0].set_xlabel('# Pathogenic Variants')
    axes[1].set_xlabel('# Benign Variants')
    for ax in axes.flat:
        ax.grid(False)
        ax.margins(0.00)
        ax.tick_params(axis='both', which='major', labelsize=8)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_tools(data, fname):
    df = data.iloc[::-1]
    set_style()

    correct_p = df['tp'].values
    incorrect_p = df['scored_p']  # [ x+y for (x,y) in zip(df['tp'].values, df['fn'].values)]
    correct_n = df['tn']
    incorrect_n = df['scored_n']  # [ x+y for (x,y) in zip(df['tn'], df['fp'])]

    w = 0.8
    fig, axes = plt.subplots(ncols=2, sharey=True)
    ind = np.arange(df.shape[0])

    axes[0].barh(ind, df['tp'], align='center', color='green', alpha=0.8, zorder=10, height=w, edgecolor='black',
                 linewidth=0.75)
    axes[0].barh(ind, df['fn'], left=correct_p, align='center', alpha=0.8, color='darkred', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)
    axes[0].barh(ind, df['mp'], left=incorrect_p, align='center', color='lightgrey', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)
    axes[1].barh(ind, df['tn'], align='center', color='green', alpha=0.8, zorder=10, height=w, edgecolor='black',
                 linewidth=0.75)
    axes[1].barh(ind, df['fp'], left=correct_n, align='center', alpha=0.8, color='darkred', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)
    axes[1].barh(ind, df['mn'], left=incorrect_n, align='center', color='lightgrey', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)

    axes[0].invert_xaxis()
    axes[0].set(yticks=ind, yticklabels=df['tool'])

    axes[0].set_xlabel('# Pathogenic Variants')
    axes[1].set_xlabel('# Benign Variants')
    for ax in axes.flat:
        ax.grid(False)
        ax.margins(0.00)
        ax.tick_params(axis='both', which='major', labelsize=8)
        for tick in ax.xaxis.get_minor_ticks():
            tick.tick1line.set_markersize(0)
            tick.tick2line.set_markersize(0)

    set_size(fig, len(df['tool']))
    fig.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


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
    for vartype, vartypefunction in filtes_var_type:

        ensure_folder_exists(os.path.join(folder, "figures", vartype))
        outdir = os.path.join(folder, "figures", vartype)
        df_v = vartypefunction(dataset).copy()
        print("\n\nLooking at {} ({} variants)".format(vartype, df_v.shape[0]))
        if df_v.shape[0] == 0:
            print("There are no {} in the variant set. Skipping this analysis.".format(vartype))
            continue

        for filtername, filterfunction in filters:
            statistics = defaultdict(list)
            df = filterfunction(df_v).copy()
            if df.shape[0] < 10:
                print("{} has not a minimum number of {} {} variants (10) to evaluate tools performance. ({})".format(
                    name, vartype, filtername, df.shape[0]))
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
