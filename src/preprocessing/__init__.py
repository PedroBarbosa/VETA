import os.path
from .utils import *
from .clinvar import *
from .vcf_processing import *
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')


def preprocess(loc, ROOT_DIR, thresholdAnalysis, clinvarStars, intronic_analysis, dataset= None, newDataset=None):
    DATASET_FOLDER = os.path.join(os.path.dirname(ROOT_DIR), "datasets")

    dict = {}
    if not dataset or dataset == "clinvar" or thresholdAnalysis:
        if thresholdAnalysis:
            logging.info("Processing clinvar data for threshold analysis..")
        else:
            logging.info("Processing clinvar data..")

        fclinvar= [filename for filename in os.listdir(os.path.join(DATASET_FOLDER, 'clinvar')) if
                   fnmatch.fnmatch(filename, "*clinvar*gz")]
        try:
            file_name = os.path.join(DATASET_FOLDER, 'clinvar', fclinvar[0])
        except IndexError:
            logging.error("Perhaps your clinvar file does not have 'clinvar' string on its name")
            exit(1)
        df_clinvar = get_clinvar_cached(file_name, loc, intronic_analysis)
        df_clinvar_sure = filter_clinvar_sure(df_clinvar)

        dict = {
            'clinvar': df_clinvar_sure,
            '1s': filter_clinvar_1_star(df_clinvar_sure),
            '2s': filter_clinvar_2_stars(df_clinvar_sure),
            '3s': filter_clinvar_3_stars(df_clinvar_sure),
            '4s': filter_clinvar_4_stars(df_clinvar_sure),
            'clinvar_l': df_clinvar,
            '1s_l': filter_clinvar_1_star(df_clinvar),
            '2s_l': filter_clinvar_2_stars(df_clinvar),
            '3s_l': filter_clinvar_3_stars(df_clinvar),
            '4s_l': filter_clinvar_4_stars(df_clinvar),
        }

    if dataset and dataset != "clinvar":
        logging.info("Prepocessing {} data".format(dataset))
        if newDataset:
            if os.path.isdir(dataset):
                dirname = dataset
                dataset_name = utils.check_dataset_arg(dataset)
            else:
                dirname = os.path.join(DATASET_FOLDER, dataset)
                dataset_name = dataset

        else:
            dirname = os.path.join(DATASET_FOLDER, dataset)
            dataset_name = dataset

        dfs = []
        benign, deleterious = check_required_files(dirname, dataset_name)
        for f in [benign, deleterious]:
            isbenign = False if "pathogenic" in f or "deleterious" in f else True
            dfs.append(get_df_ready(f, False, isbenign, loc, intronic_analysis))

        df = pd.concat(dfs)
        df = vcf_cleaning(df)
        dict[dataset_name] = df

    return dict
