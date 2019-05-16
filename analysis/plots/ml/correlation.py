import os

import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

plt.switch_backend('agg')

def generate_ml_feature_correlation(df, name, threshold_list, folder):
    """Generates Pearson's feature correlation"""
    plt.figure(figsize=(18, 12))
    sns.heatmap(df[[tool[0] for tool in threshold_list]].corr(), cmap="YlGnBu",annot=False, fmt=".2f")
    plt.tight_layout()
    plt.savefig(os.path.join(folder, 'figures','pearson_feat_corr_{}.pdf'.format(name)))
    plt.close()