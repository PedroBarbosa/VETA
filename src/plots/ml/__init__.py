from .plot_feature_correlation import generate_ml_feature_correlation
from predictions.ml_feature_ranking import generate_feature_ranking
from predictions.ml_classifiers import generate_classifiers_analysis
from predictions.ml_prepare_dataset import prepare_dataset, undersample
from .plot_tree_builder import generate_tree_plot
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

def generate_ml_analysis(df, filters, threshold_list, name, folder):
    
    for filtern, filterf in filters:
        df_ = filterf(df).copy()
        X, y = prepare_dataset(df_, threshold_list)
        X, y = undersample(X, y)
        
        if len(X) < 10:
            logging.info("{} {} has less than 10 samples. Skipping.".format(name, filtern))
            continue
        
        featsel, new_thresholds = generate_feature_ranking(df_.columns, X, y, threshold_list, filtern, name, folder)
        X_small = featsel.transform(X)
        generate_classifiers_analysis(X, y, threshold_list, filtern, name, folder)
        generate_classifiers_analysis(X_small, y, new_thresholds, filtern, name, folder, extra="less_features")
        generate_tree_plot(X_small, y, new_thresholds, filtern, name, folder)