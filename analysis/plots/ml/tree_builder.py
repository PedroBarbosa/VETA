import os.path

from sklearn.datasets import *
from sklearn import tree
from dtreeviz.trees import *

def generate_tree_plot(X, y, threshold_list, filtern, name, folder):
    classifier = tree.DecisionTreeClassifier(
        criterion = "entropy",
#        min_samples_leaf = 5,
        presort = True,
        max_depth = 4)
    classifier.fit(X, y)

    viz = dtreeviz(classifier,
                   X,
                   y,
                   target_name = 'pathogenicity',
                   feature_names = [ x[0] for x in threshold_list],
                   class_names = ['benign', 'pathogenic', 'unknown'])
              
    viz.save(os.path.join(folder, "figures", "tree_{}_{}.svg".format(name, filtern)))