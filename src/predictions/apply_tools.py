import numpy as np
import pandas as pd
import logging
import sys
logging.basicConfig(stream=sys.stdout, level=logging.INFO,  format='%(asctime)s %(message)s')

def apply_tool_predictions(df, thresholds):
    for tool, direction, threshold, *args in thresholds:
        if tool not in df.columns:
            logging.info(tool + " scores were not found. Probably they are not present in the input VCF.")
            continue

        if direction == ">":
            classification = lambda x: pd.isna(x) and np.nan or x > threshold
        else:
            classification = lambda x: pd.isna(x) and np.nan or x < threshold
        
        df[tool + '_prediction'] = df[tool].apply(classification)
    return df
