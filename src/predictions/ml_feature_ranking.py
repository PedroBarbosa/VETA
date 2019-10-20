import os

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')
sns.set(style="white")
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.model_selection import StratifiedKFold
from predictions.ml_cost_aware_rfecv import CostAwareRFECV
from yellowbrick.model_selection import RFECV
from plots.ml.plot_feature_correlation import plot_feature_importance


def information_gain(X, y, threshold_list, filtern, folder):
    ig = mutual_info_classif(X, y, discrete_features=False)
    indices = np.argsort(ig)[::-1]
    sorted_tools = [ threshold_list[i][0] for i in indices ]
    plot_feature_importance(sorted_tools, X, ig, indices, [0 for i in indices], filtern, 'information_gain', folder)


def feature_importance(X, y, threshold_list, filtern, folder):
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

    plot_feature_importance(sorted_tools, X, importances, indices, std, filtern, 'feature_importance', folder)


def recursive_feature_elimination(X, y, threshold_list, filtern, folder):
    """This method is wrongly printing the order of the ranked features. Need to fix it.
     Also it should detect when a smaller number of features is enough to achieve almost as much as max.
     Right now, picks the max value"""
    
    rf = RandomForestClassifier(n_estimators=10, criterion='entropy')
    image_path = os.path.join(folder,  "rfecv_{}.pdf".format(filtern))
    selection_path = os.path.join(folder,  "rfecv_{}.tex".format(filtern))

    fig = plt.figure()
    rfecv = RFECV(rf, cv=StratifiedKFold(5), scoring='f1_weighted', step=1, error=0.15)
    #rfecv = CostAwareRFECV(rf, cv=StratifiedKFold(5), scoring='f1_weighted', step=1, error=0.15)
    rfecv.fit(X, y)
    rfecv.poof(outpath=image_path)
    tools = [x[0] for x in threshold_list]
    selected_idx = [i for i, x in enumerate(rfecv.support_) if x]

    selected = [tools[i] for i in selected_idx]
    with open(selection_path, "w") as f:
        f.write(", ".join(selected))
    plt.close()

    return rfecv, [threshold_list[i] for i in selected_idx]


def generate_feature_ranking(columns, X, y, threshold_list, filtern, folder):
    feature_importance(X, y, threshold_list, filtern, folder)
    information_gain(X, y, threshold_list, filtern, folder)
    adapter, new_features = recursive_feature_elimination(X, y, threshold_list, filtern, folder)
    return adapter, new_features
