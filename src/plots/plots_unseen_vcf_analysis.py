import matplotlib.pyplot as plt
import os
import seaborn as sns
import matplotlib.patches as mpatches

def plot_area(df, outdir):
    fig, ax = plt.subplots\
        (figsize=(12,12))
    map_colors = {"unpredictable": "grey",
              "is_benign": "darkblue",
              "is_pathogenic": "darkred"
              }

    ax = df.plot.area(color=[map_colors.get(x) for x in df.columns], alpha=0.85)
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    ax.set_ylabel('Fraction of tools')
    ax.set_xlabel('Variants')
    plt.xticks([])
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, 'areaplot.pdf'))
    plt.close()

def plot_heatmap(df, outdir, displayAnnot):
    fig, ax = plt.subplots(1, 2, figsize=(5,7))
    #cbaxes = fig.add_axes([0.2, 0, 0, 1])
    #sns.heatmap(df[["is_pathogenic"]], ax=ax[0], cmap="PiYG_r", cbar_ax=cbaxes)#cbar_kws = dict(cax=cbaxes))
    sns.heatmap(df[["is_pathogenic"]], ax=ax[0], vmax=1, vmin=0, cmap="OrRd", yticklabels=displayAnnot)
    sns.heatmap(df[["unpredictable"]], ax=ax[1], cbar_kws={'label': '%'}, vmax=100, vmin=0,  cmap="OrRd", yticklabels= False)
    ax[0].set_ylabel('')
    ax[1].set_ylabel('')
    plt.tight_layout()
    if displayAnnot:
        plt.savefig(os.path.join(outdir, 'heatmap_topvariants.pdf'))
    else:
        plt.savefig(os.path.join(outdir, 'heatmap_all.pdf'))
    plt.close()


def plot_heatmap_toptools(df,filters,outdir):
    for filtername, filterfunction in filters:
        df_f = filterfunction(df).drop(["location"], axis=1).copy()
        if df_f.shape[0] > 0:
            plt.clf()
            f, ax = plt.subplots(figsize=(0.5 * len(df_f.columns), 6))
            colors = ["slateblue", "w", "indianred"]
            sns.heatmap(df_f, cmap=colors, ax=ax, cbar=False)
            if df_f.shape[0] > 40:
                plt.yticks([])
                ax.set_ylabel("{} variants".format(df_f.shape[0]))
            else:
                ax.set_ylabel('')
            box = ax.get_position()
            ax.set_position([box.x0, box.y0, box.width * 0.7, box.height])
            legend_ax = f.add_axes([.7, .5, 1, .1])
            legend_ax.axis('off')

            patches = [mpatches.Patch(facecolor=c, edgecolor=c) for c in colors]
            legend_ax.legend(patches, ["Benign", "Unpredictable", "Pathogenic"], handlelength=0.8, loc='lower left')
            plt.savefig(os.path.join(outdir,"top_tools_{}.pdf".format(filtername)),bbox_inches='tight', pad_inches = 0)
            plt.close()

def plot_tool_score_distribution(df, tool, threshold_list, outdir):
   # sns.kdeplot(df[tool], shade=True, cut=0)
   # sns.rugplot(df[tool])
    df[tool].dropna(inplace=True)
    sns.distplot(df[tool], hist=True, kde=False,
                bins=30, color='slateblue',
                hist_kws={'edgecolor': 'black'}
                )

    plt.axvline([x[2] for x in threshold_list if x[0] == tool], color='r', linestyle="--")
    plt.savefig(os.path.join(outdir,"score_dist_{}.pdf".format(tool)))
    plt.close()
