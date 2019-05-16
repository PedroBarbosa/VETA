import os.path
from .utils import *
from .clinvar import *
from .vep import *
DATASET_FOLDER = "../datasets/"


def preprocess(loc, dataset=None):

    fclinvar= [filename for filename in os.listdir(os.path.join(DATASET_FOLDER,'clinvar')) if fnmatch.fnmatch(filename,"clinvar*gz")]
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

    if dataset:
        dfs=[]
        benign, deleterious = check_required_files(os.path.join(DATASET_FOLDER, dataset), dataset)
        for f in [benign,deleterious]:
            isbenign= False if "pathogenic" in f or "deleterious" in f else True
            dfs.append(get_df_ready(f, False, isbenign, loc))

        df = pd.concat(dfs)
        df = vep_cleaning(df)
        dict[dataset] = df

    return dict