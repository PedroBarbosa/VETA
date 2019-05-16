import copy
import os.path

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
plt.switch_backend('agg')
sns.set(style="white")

from sklearn.metrics import confusion_matrix, accuracy_score
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline


from sklearn import tree
from sklearn import dummy
from sklearn import svm
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier, AdaBoostClassifier
from sklearn.linear_model import SGDClassifier
from sklearn.linear_model.logistic import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.neighbors.nearest_centroid import NearestCentroid

from lightning.classification import CDClassifier
from pyearth import Earth
from gplearn.genetic import SymbolicTransformer
from xgboost import XGBClassifier
from bleedml.classifiers import CascadeForest

from .single_feature_classifier import SingleFeatureClassifier
from ..utils import *

classifiers = [
#    ('Random', dummy.DummyClassifier()),    
    ('Naive Bayes', GaussianNB()),
    ('QDA', QuadraticDiscriminantAnalysis()),    
    ('KNN', NearestCentroid()),
    ('Decision Tree', tree.DecisionTreeClassifier()),
    ('Extra Trees', ExtraTreesClassifier()),
    ('Random Forests', RandomForestClassifier(n_estimators=100)),
    ('Ada Boost', AdaBoostClassifier()),
    ('SVM OVO', svm.SVC(decision_function_shape='ovo', gamma='scale')),
    ('MLP', MLPClassifier(solver='lbfgs', alpha=1e-5)),
    ('Lightining', CDClassifier(penalty="l1/l2",
                   loss="squared_hinge",
                   multiclass=False,
                   max_iter=20)),
    ('Earth', Pipeline([('earth', Earth(max_degree=3, penalty=1.5)),
                             ('logistic', LogisticRegression(solver='lbfgs'))])),
    ('GP', Pipeline([('gp', SymbolicTransformer()),
                             ('logistic', LogisticRegression(solver='lbfgs'))])),
    ('XGBoost', XGBClassifier(n_jobs = -1)),
    ('Cascade Forests', CascadeForest())

]

def expand_classifier_list(threshold_list):
    classifiers_total = copy.deepcopy(classifiers)
    i = 0
    for t, a, b, c, m in threshold_list:
        single_cls = SingleFeatureClassifier(t, a, b, i)
        classifiers_total.append((t, single_cls))
        i += 1
    return classifiers_total

def make_metric(yreal, ypred):
    return {
        'accuracy': accuracy_score(yreal, ypred),
        'tp': len([ 1 for (r,p) in zip(yreal, ypred) if r == 1 and p == r ]),
        'tn': len([ 1 for (r,p) in zip(yreal, ypred) if r == 0 and p == r ]),
        'fp': len([ 1 for (r,p) in zip(yreal, ypred) if r == 0 and p == 1 ]),
        'fn': len([ 1 for (r,p) in zip(yreal, ypred) if r == 1 and p == 0 ]),
        'mp': len([ 1 for (r,p) in zip(yreal, ypred) if r == 1 and p not in [0,1] ]),
        'mn': len([ 1 for (r,p) in zip(yreal, ypred) if r == 0 and p not in [0,1] ])
        }

def perform_cross_validation(X, y, name, classifier):
    gkf = StratifiedKFold(10)
    
    for train, test in gkf.split(X, y):
        classifier.fit(X[train], y[train])
        yreal = y[test]
        ypred = classifier.predict(X[test])
        score = make_metric(yreal, ypred)
        yield score
    

def plot_classifiers(out_metrics, filtern, fname, folder, extra=None):
    fig, ax = plt.subplots(figsize=(15,8))
        
    classifier_names = [ k for k in out_metrics.keys() if np.average([p['accuracy'] for p in out_metrics[k]]) > 0 ]
    
    ind = np.arange(len(classifier_names))
    
    correct_p = [ np.average([p['tp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    incorrect_p = [ np.average([p['fn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    correct_n = [ np.average([p['tn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    incorrect_n = [ np.average([p['fp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    
    missing_p = [ np.average([p['mp'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    missing_n = [ np.average([p['mn'] for p in out_metrics[cl_name]]) for cl_name in classifier_names]
    
    w = 0.8
    
    fig, axes = plt.subplots(ncols=2, sharey=True)
    axes[0].barh(ind, correct_p, align='center', color='green', zorder=10, height=w, linewidth=0)
    axes[0].barh(ind, incorrect_p, left=correct_p, align='center', color='darkred', zorder=10, height=w, linewidth=0)
    axes[0].barh(ind, missing_p, left=[ x+y for (x,y) in zip(incorrect_p, correct_p)], align='center', color='lightgrey', zorder=10, height=w, linewidth=0)
    axes[1].barh(ind, correct_n, align='center', color='green', zorder=10, height=w, linewidth=0)
    axes[1].barh(ind,incorrect_n, left=correct_n, align='center', color='darkred', zorder=10, height=w, linewidth=0)
    axes[1].barh(ind, missing_n, left=[ x+y for (x,y) in zip(incorrect_n, correct_n)], align='center', color='lightgrey', zorder=10, height=w, linewidth=0)

    axes[0].invert_xaxis()
    axes[0].set(yticks=ind, yticklabels=["{} ({:.0f}%)".format(c, np.average([p['accuracy'] for p in out_metrics[c]]) * 100) for c in classifier_names])
        
    axes[0].set_xlabel('Pathogenic Variants')
    axes[1].set_xlabel('Benign Variants')

    set_size(fig, len(classifier_names))    
    fig.tight_layout()
    if extra:
        extra = "_" + extra
    else:
        extra = ""
    plt.savefig(os.path.join(folder, "figures", "classifiers_comparison_{}_{}{}.pdf".format(fname, filtern, extra)))

def generate_classifiers_analysis(X, y, threshold_list, filtern, fname, folder, extra=None):
    if extra:
        classifiers = expand_classifier_list(threshold_list)
    else:
        classifiers = expand_classifier_list(threshold_list)
    
    out_metrics = {}
    for name, classifier in classifiers:
        out_metrics[name] = list(perform_cross_validation(X, y, name, classifier))
        
    plot_classifiers(out_metrics, filtern, fname, folder, extra)