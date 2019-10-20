import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')


def generate_ml_feature_correlation(df, threshold_list, outdir):
    """Generates Pearson's feature correlation"""
    plt.figure(figsize=(18, 12))
    sns.heatmap(df[[tool[0] for tool in threshold_list if tool[0] in df.columns]].corr(), cmap="YlGnBu",annot=False, fmt=".2f")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'pearson_feat_corr.pdf'))
    plt.close()


def plot_feature_importance(sorted_tools, X, importances, indices, std, filtern, method, folder):
    col = sns.color_palette("Blues_r",10)
    while len(col) < len(sorted_tools):
        col.append(matplotlib.colors.to_rgba('grey'))
    fig = plt.figure()

    if method == 'information_gain':
        p = plt.barh(range(X.shape[1]), importances[indices],
                      color=col, align="center")
    elif method == 'feature_importance':
        p = plt.barh(range(X.shape[1]), importances[indices], yerr=std[indices], color=col,  align="center")
        plt.xlabel('Feature Importance')

    #plt.legend(fig, sorted_tools, loc="upper right")
    plt.xlim(left=0)
    plt.ylim([-1, len(sorted_tools)])
    plt.yticks(range(0,len(sorted_tools)), sorted_tools)
    fig.tight_layout()
    plt.savefig(os.path.join(folder,  "feature_ml_{}_{}.pdf".format(filtern, method)))
    plt.close()
