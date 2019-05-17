import os.path
from .utils import *
from .clinvar import *
from .vep import *
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

def preprocess(loc, ROOT_DIR, thresholdAnalysis, dataset=None, newDataset=None):
    DATASET_FOLDER = os.path.join(os.path.dirname(ROOT_DIR), "datasets")

    dict = {}
    if not dataset or dataset == "clinvar" or thresholdAnalysis:
        if thresholdAnalysis:
            logging.info("Processing clinvar data for threshold analysis..")
        else:
            logging.info("Processing clinvar data..")

        fclinvar= [filename for filename in os.listdir(os.path.join(DATASET_FOLDER, 'clinvar')) if
                   fnmatch.fnmatch(filename, "clinvar*gz")]
        file_name = os.path.join(DATASET_FOLDER, 'clinvar', fclinvar[0])
        df_clinvar = get_clinvar_cached(file_name,loc)
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
            dirname = dataset
            if dataset.rsplit('/', 1)[-1]:
                dataset_name = dataset.rsplit('/', 1)[-1]
            else:
                dataset_name = dataset.rsplit('/', 1)[-2]

        else:
            dirname = os.path.join(DATASET_FOLDER, dataset)
            dataset_name = dataset

        dfs = []
        benign, deleterious = check_required_files(dirname, dataset_name)
        for f in [benign,deleterious]:
            isbenign = False if "pathogenic" in f or "deleterious" in f else True
            dfs.append(get_df_ready(f, False, isbenign, loc))

        df = pd.concat(dfs)
        df = vep_cleaning(df)
        dict[dataset] = df

    return dict
