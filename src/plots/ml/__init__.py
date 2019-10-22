from .plot_feature_correlation import generate_ml_feature_correlation
from predictions.ml_feature_ranking import generate_feature_ranking
from predictions.ml_classifiers import generate_classifiers_analysis
from predictions.ml_prepare_dataset import prepare_dataset, undersample
from .plot_tree_builder import generate_tree_plot
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


def generate_ml_analysis(df, filters, threshold_list, name, folder):

    for filtern, filterf in [filt for filt in filters if filt[0] in ["all", "coding", "splicesite"]]:
        df_ = filterf(df).copy()
        X, y, threshold_list_present = prepare_dataset(df_, threshold_list)
        X, y = undersample(X, y)
        
        if len(X) < 10:
            logging.info("{} {} has less than 10 samples. Skipping.".format(name, filtern))
            continue

        logging.info("Looking at {} variants.".format(filtern))
        logging.info("Generating features correlation")
        generate_ml_feature_correlation(df_, threshold_list, filtern, folder)
        logging.info("Generating feature ranking.")
        featsel, tools_subset = generate_feature_ranking(df_.columns, X, y, threshold_list_present, filtern, folder)
        X_small = featsel.transform(X)
        logging.info("Applying classifiers to data.")
        generate_classifiers_analysis(X, y, threshold_list_present, filtern, folder)
        logging.info("Applying classifiers to data with less features..")
        logging.info("Tools selected: {}".format([t[0] for t in tools_subset]))
        generate_classifiers_analysis(X_small, y, tools_subset, filtern, folder, extra="less_features")
        generate_tree_plot(X_small, y, tools_subset, filtern, name, folder)
