import os.path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys
import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')
plt.switch_backend('agg')
sns.set(style="white")
GT_LABEL = '(*) Ground Truth'

def plot_heatmap_unlabelled(df):
    f, ax = plt.subplots(figsize=(0.5 * len(df.columns), 6))
    cmap = sns.diverging_palette(10, 220, as_cmap=True)
    sns.heatmap(df, cmap=cmap,  linewidths=0, cbar=False)
    plt.show()
    
def prepare_dataset_for_heatmap(df, thresholds):
    df['blank'] = pd.Series(np.nan, index = np.arange(df.shape[0]))
    df = df[['class'] + [col for col in df.columns if '_prediction' in col]].copy()
    df.columns = [GT_LABEL] + [col.replace("_prediction", "") for col in df.columns if '_prediction' in col]
    f = lambda x: pd.isna(x) and np.nan or (1-int(x))
    for col in df.columns:
        df[col] = df[col].apply(f)
        if df[col].isnull().all():
            del df[col]
    df = df.reset_index(drop=True)
    df = df.sort_values(by=GT_LABEL)
    return df

def hide_ylabel():
    fr = plt.gca()
    fr.axes.yaxis.set_ticklabels([])

def plot_heatmap(df, filtername, thresholds, fname):

    df = prepare_dataset_for_heatmap(df, thresholds)
    plt.clf()
    f, ax = plt.subplots(figsize=(0.5 * len(df.columns), 6))
    cmap="PiYG"
    cmap = sns.diverging_palette(10, 220, as_cmap=True)
    sns.heatmap(df, cmap=cmap, linewidths=0, cbar=False, mask=df.isnull())
    hide_ylabel()
    plt.savefig(fname, bbox_inches='tight')
    plt.close()


def generate_heatmap(df_original, filters_var_type, filters, thresholds, name, folder):
    logging.info("-----------------------------")
    logging.info("Generating performance heatmaps. If dataset is large, this step will take long time.")
    logging.info("-----------------------------")
    for vartype, vartypefunction in filters_var_type:
        if not os.path.exists(os.path.join(folder, "tools_benchmark", vartype)):
            os.mkdir(os.path.join(folder, "tools_benchmark", vartype))
        outdir = os.path.join(folder, "tools_benchmark", vartype)
        df_v = vartypefunction(df_original).copy()
        for filtername, filterfunction in filters:
            df = filterfunction(df_v).copy()
            if df.shape[0] > 10:
                plot_heatmap(df, filtername, thresholds, os.path.join(outdir, "heatmap_{}_{}.pdf".format(name,
                                                                                            filtername)))
