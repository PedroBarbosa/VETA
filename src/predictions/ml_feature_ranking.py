import os

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')
sns.set(style="white")
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import StratifiedKFold
from predictions.ml_cost_aware_rfecv import CostAwareRFECV


def plot_feature_importance(sorted_tools, X, importances, indices, std, filtern, name, method, folder):
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
    plt.savefig(os.path.join(folder, "figures", "feature_ml_{}_{}_{}.pdf".format(name, filtern, method)))
    plt.close()
    
    
def information_gain(X, y, threshold_list, filtern, name, folder):
    ig = mutual_info_classif(X,y,discrete_features=False)
    indices = np.argsort(ig)[::-1]
    sorted_tools = [ threshold_list[i][0] for i in indices ]
    plot_feature_importance(sorted_tools, X, ig, indices, [0 for i in indices], filtern, name, 'information_gain', folder)

def feature_importance(X, y, threshold_list, filtern, name, folder):
    forest = RandomForestClassifier(n_estimators=100, criterion='entropy')
    forest.fit(X, y)
    importances = forest.feature_importances_
    tools=[x[0] for x in threshold_list]

    std = np.std([tree.feature_importances_ for tree in forest.estimators_],
                 axis=0)
    indices = np.argsort(importances)[::-1]

    sorted_tools=[]
    for f in range(X.shape[1]):
        if importances[indices[f]] > 0:
            sorted_tools.append(tools[indices[f]])

    plot_feature_importance(sorted_tools, X, importances, indices, std, filtern, name, 'feature_importance', folder)


def recursive_feature_elimination(X, y, threshold_list, filtern, name, folder):
    """This method is wrongly printing the order of the ranked features. Need to fix it. Also it should detect when a smaller number of features is enough to achieve almost as much as max.
     Right now, picks the max value"""
    
    rf = RandomForestClassifier(n_estimators=10, criterion='entropy')

    image_path = os.path.join(folder, "figures", "rfecv_{}_{}.pdf".format(name,filtern))
    selection_path = os.path.join(folder, "figures", "rfecv_{}_{}.tex".format(name,filtern))

    fig = plt.figure()
    rfecv = CostAwareRFECV(rf, cv=StratifiedKFold(5), scoring='f1_weighted', step=1, error=0.15)
    rfecv.fit(X, y)
    rfecv.poof(outpath=image_path)
    tools = [x[0] for x in threshold_list]
    selected_idx=[i for i,x in enumerate(rfecv.support_) if x]
    selected=[tools[i] for i in selected_idx]
    with open(selection_path, "w") as f:
        f.write(", ".join(selected))
    plt.close()

    return rfecv, [threshold_list[i] for i in selected_idx]

def generate_feature_ranking(columns, X, y, threshold_list, filtern, name, folder):
    feature_importance(X, y, threshold_list, filtern, name, folder)
    information_gain(X, y, threshold_list, filtern, name, folder)
    adapter, new_features = recursive_feature_elimination(X, y, threshold_list, filtern, name, folder)
    return adapter, new_features
