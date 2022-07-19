import pandas as pd
import numpy as np

from ukbb_data_scripts import constants


def get_features_df(df, features):
    out_df = df.copy()
    # keep necessary fields
    out_df = out_df[["var_id", "predictand", "highest_AF"]]
    for feature in features:
        # simple features
        if feature in set(df.columns.tolist()):
            if feature in constants.all_cat_cols:
                cat_dummies = pd.get_dummies(df[feature]).values
                for i in range(cat_dummies.shape[1]):
                    out_df[feature + f"_{i}"] = cat_dummies[:, i]
            elif feature in constants.all_cont_cols:
                out_df[feature] = df[feature]
                out_df[feature + "_na"] = df[feature + "_na"]
        else:
            # this clause corresponds to either a power or log of the feature
            orig_feature, operand = feature.rsplit("_", 1)
            out_df[orig_feature + "_na"] = df[orig_feature + "_na"]
            if operand == "**2":
                out_df[feature] = (df[orig_feature]**2).values
            elif operand == "log":
                out_df[feature] = np.log(np.abs(df[orig_feature]) + 1).values
            else:
                print(f"feature operand not recognized: {operand}")
    return out_df
