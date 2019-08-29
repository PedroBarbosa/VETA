import numpy as np
import seaborn as sns
sns.set(style="white")
from preprocessing.utils import ratio
import os
from collections import defaultdict
from osutils import ensure_folder_exists
from matplotlib import lines
from plots.plots_utils import *
import matplotlib.pyplot as plt
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


def plot_unscored(data, fname):
    ax = sns.barplot(x="fraction_nan", y="tool", color="skyblue", data=data)
    ax.set(xlabel='Fraction of unscored variants', ylabel='')
    plt.tight_layout()
    plt.savefig(fname + ".pdf")
    plt.close()


def plot_precision_recall(data, threshold_list, fname):
    df_precision_recall = pd.DataFrame(columns=['tool','threshold', 'precision', 'recall'])
    for tool, direction, recommended_threshold, *args in threshold_list:
        try:
            df_ = data.loc[pd.notnull(data[tool]), ].copy()
        except KeyError:
            continue

        max_thr = df_[tool].max()
        min_thr = df_[tool].min()

        if pd.isnull(max_thr) or pd.isnull(min_thr) or max_thr == min_thr:
            print("Something strange in max/min thresholds {} {} {}".format(tool, max_thr, min_thr))
            continue

        step = (max_thr - min_thr) / float(100)
        threshold_range = np.arange(min_thr, max_thr, step)

        toappend=[]
        for threshold in threshold_range:
            if direction == ">":
                classification_f = lambda x: x == np.nan and np.nan or x > threshold
            else:
                classification_f = lambda x: x == np.nan and np.nan or x < threshold

            classification = df_[tool].map(classification_f)
            # df_ = df.copy()

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
            f1 = ratio(2.0 * (precision * sensitivity), (sensitivity + precision))


            toappend.append([tool, threshold, precision, sensitivity])
        pd.concat([df_precision_recall, pd.Series(toappend)], ignore_index=True)

    print(df_precision_recall.head)


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

    ax.grid(axis='x', linestyle='dashed')
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

    axes[0].barh(ind, df['tp'], align='center', color='darkblue', alpha=0.7, zorder=10, height=w, edgecolor='black',
                 linewidth=0.75)
    axes[0].barh(ind, df['fn'], left=correct_p, align='center', alpha=0.8, color='darkred', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)
    axes[0].barh(ind, df['mp'], left=incorrect_p, align='center', color='lightgrey', zorder=10, height=w,
                 edgecolor='black', linewidth=0.75)
    axes[1].barh(ind, df['tn'], align='center', color='darkblue', alpha=0.7, zorder=10, height=w, edgecolor='black',
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
