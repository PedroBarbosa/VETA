import numpy as np
import pandas as pd

def apply_tool_predictions(df, thresholds):
    for tool, direction, threshold, color, marker in thresholds:
        if tool not in df.columns:
            print("Problem with ", tool)
            return
        if direction == ">":
            classification = lambda x: pd.isna(x) and np.nan or x > threshold
        else:
            classification = lambda x: pd.isna(x) and np.nan or x < threshold
        
        df[tool + '_prediction'] = df[tool].apply(classification)
    return df
