import pandas as pd
import seaborn as sns
from dtreeviz.trees import *
from sklearn.feature_selection import _rfe
from predictions.utils_classifiers import *


def plot_feature_correlation(df: pd.DataFrame, 
                             location: str, 
                             outdir: str,
                             drop_AF: bool = True):
    """
    Generates Pearson's feature correlation

    :param pd.DataFrame df: Input dataframe
    :param str location: Locations of the variants
    :param str outdir: Output directory
    """
    cmap = sns.diverging_palette(220, 10, as_cmap=True)

    if "gnomADg_AF" in df.columns and drop_AF:
        df.drop('gnomADg_AF', axis=1, inplace=True)

    if 'intron_offset' in df.columns:
        corr = df.drop(['label', 'intron_offset'], axis=1).corr()
    else:
        corr = df.drop('label', axis=1).corr()

    # Generate a mask for the upper triangle
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True

    annot=False if len(corr.columns) > 15 else True
    f, ax = plt.subplots(figsize=(11, 9))
    ax = sns.heatmap(corr,
                     mask=mask,
                     cmap=cmap,
                     vmax=1,
                     vmin=-1,
                     annot=annot, 
                     square=True,
                     linewidths=.5,
                     linecolor='k',
                     cbar_kws={'shrink': 0.5},
                     ax=ax)

    ax.axhline(y=0, color='k', linewidth=5)
    ax.axhline(y=corr.shape[1], color='k', linewidth=5)
    ax.axvline(x=0, color='k', linewidth=5)
    ax.axvline(x=corr.shape[0], color='k', linewidth=10)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=70)
    plt.tight_layout()
    outdir = os.path.join(outdir, 'feature_correlation')
    os.makedirs(outdir, exist_ok=True)
    plt.savefig(os.path.join(outdir, 'pearson_correlation_{}.pdf'.format(location)))
    plt.close()


def plot_feature_importance(importance_scores: dict,
                            method: str,
                            location: str,
                            out_dir: str,
                            std: np.array = None):
    """
    Plot feature importance scores

    :param dict importance_scores: Feature names with corresponding
        importance scores
    :param str method: Whether scores are from random forests
    `feature_importance` attribute or from information gain analysis
    :param str location: Location of the variants
    :param str out_dir: Output directory
    :param np.array std: Standard deviation for the importance scores.
    """
    fig = plt.figure()

    # Plot just top 20
    if len(importance_scores) > 20:
        importance_scores = {k: v for (k, v) in [x for x in importance_scores.items()][-20:]}

        if std is not None:
            std = std[-20:]

    if method == 'feature_importance':

        p = plt.barh(range(len(importance_scores)), importance_scores.values(),
                     xerr=std,
                     error_kw={'markeredgewidth': 10},
                     capsize=1,
                     color='darkgrey',
                     edgecolor='k',
                     linewidth=1)
        plt.xlabel('Feature Importance')

    elif method == 'information_gain':
        out_dir = os.path.join(out_dir, "mutual_information_gain")
        os.makedirs(out_dir, exist_ok=True)
        p = plt.barh(range(len(importance_scores)), importance_scores.values(),
                     color='darkgrey',
                     edgecolor='k',
                     linewidth=1)
        plt.xlabel('Information Gain')

    plt.xlim(left=0)
    plt.yticks(range(0, len(importance_scores)), importance_scores.keys())
    plt.ylim([-1, len(importance_scores)])
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, "{}_{}.pdf".format(method, location)))
    plt.close()


def plot_rfecv(selector: _rfe.RFECV,
               min_features: int,
               feature_names: list,
               clf_name: str,
               loc: str,
               out_dir: str):
    """
    Plot cross-validation performance after recursive feature elimination

    :param _rfe.RFECV selector: Feature selector
    :param int min_features: Minimum number of features used
    :param list feature_names: Name of the features to plot
    :param str clf_name: Name of the classifier fitted
    :param str loc: Location of the variants
    :param str out_dir: Output directory
    :return:
    """
    selected_idx = [i for i, x in enumerate(selector.support_) if x]
    selected_tools = [feature_names[i] for i in selected_idx]

    fig = plt.figure()

    plt.plot(range(min_features,
                   len(selector.cv_results_['mean_test_score']) + min_features),
             selector.cv_results_['mean_test_score'])

    plt.axvline(selector.n_features_,
                color='red',
                linestyle='--',
                linewidth=0.5)

    plt.gcf().text(0.75, 0.5,
                   'Selected features:\n{}'.format('\n'.join(selected_tools)),
                   fontsize=8)

    plt.xlabel("Number of features selected")
    plt.ylabel("Cross Validation score")
    plt.title("Number of features for best score: {}".format(selector.n_features_))
    plt.ylim(0.5, 1)
    plt.subplots_adjust(left=0.3)
    fig.tight_layout()
    plt.savefig(os.path.join(out_dir, "rfecv_{}_{}.pdf".format(loc, clf_name)))
    plt.close()


def plot_decision_tree(dt_pipeline: Pipeline,
                       x: np.array, y: np.array,
                       feature_names: list,
                       loc: str,
                       out_dir: str):
    """
    Generate decision tree decisions with dtreeviz packages
    :param dt_pipeline: Fitted estimator
    :param np.array x: Instances to run through the model
    :param np.array y: Labels
    :param list feature_names: Names of the features
    :param str loc: Location
    :param str out_dir: Output directory
    """
    _clf = dt_pipeline.named_steps['dt']
    out_dir = os.path.join(out_dir, "decision_tree")
    os.makedirs(out_dir, exist_ok=True)
    for bool_ in [False, True]:
        viz = dtreeviz(_clf,
                       x,
                       y,
                       target_name='Pathogenicity',
                       feature_names=feature_names,
                       fancy=bool_,
                       class_names=['Benign', 'Pathogenic', 'unknown'])

        out_str = "tree_{}_fancy.svg".format(loc) if bool_ else "tree_{}.svg".format(loc)
        viz.save(os.path.join(out_dir, out_str))
