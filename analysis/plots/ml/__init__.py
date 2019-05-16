import numpy as np

from .correlation import generate_ml_feature_correlation
from .feature_ranking import generate_feature_ranking
from .classifiers import generate_classifiers_analysis
from .prepare_dataset import prepare_dataset, undersample
from .tree_builder import generate_tree_plot

def generate_ml_analysis(df, filters, threshold_list, name, folder):
    
    for filtern, filterf in filters:
        df_ = filterf(df).copy()
        X, y = prepare_dataset(df_, threshold_list)
        X, y = undersample(X, y)
        
        if len(X) < 10:
            print("{} {} has less than 10 samples. Skipping.".format(name, filtern))
            continue
        
        featsel, new_thresholds = generate_feature_ranking(df_.columns, X, y, threshold_list, filtern, name, folder)
        
        X_small = featsel.transform(X)
        
        generate_classifiers_analysis(X, y, threshold_list, filtern, name, folder)
        
        generate_classifiers_analysis(X_small, y, new_thresholds, filtern, name, folder, extra="less_features")
        
        generate_tree_plot(X_small, y, new_thresholds, filtern, name, folder)