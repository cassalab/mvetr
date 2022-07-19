import pandas as pd
import numpy as np

from ukbb_data_scripts import constants
from model_helpers import LogitRegression
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import global_constants

import itertools

from multiprocessing import Pool

all_vars_phenotypes_df = pd.read_pickle(global_constants.top_regression_vars_phenotypes_path)

lr_results_df = pd.read_pickle("regression_results_processed.pkl")

top_lr_results = lr_results_df.sort_values("mean", ascending=False).drop_duplicates("phenotype/gene")[["phenotype/gene", "feature", "mean", "std", "maf"]]

def process_phenotype_gene(phenotype_gene):
    train_size = 0.8
    n_seeds = 100
    phenotype = phenotype_gene.split("/")[0]
    gene = phenotype_gene.split("/")[1]
    results = []
    vars_df = all_vars_phenotypes_df[(all_vars_phenotypes_df.phenotype == phenotype) & (all_vars_phenotypes_df.gene == gene)]
    phenotype_continuous = vars_df.phenotype_continuous.iloc[0]
    # start with top result feature
    top_feature = top_lr_results[top_lr_results["phenotype/gene"] == phenotype_gene].feature.iloc[0]
    for maf in [0.1, 0.01, 0.001, 0.0001]:
        cur_vars_df = vars_df[vars_df.highest_AF < maf]
        cur_unique_vars = cur_vars_df.var_id.unique()
        cur_features_df = cur_vars_df[["var_id", "predictand"]]
        cur_best_mean_r2 = -99999
        cur_best_median_r2 = 0
        cur_best_std_r2 = 0
        cur_best_features = []

        for new_feature in ([top_feature] + constants.all_cat_cols + constants.all_cont_cols):
            features_to_add = []
            if new_feature in constants.all_cat_cols:
                features_to_add.append((new_feature, pd.get_dummies(cur_vars_df[new_feature]).values))
            else:
                features_to_add.append((new_feature, cur_vars_df[new_feature].values))
                features_to_add.append((new_feature + "_na", cur_vars_df[new_feature + "_na"].values))
                features_to_add.append((new_feature + "_**2", (cur_vars_df[new_feature]**2).values))
                features_to_add.append((new_feature + "_log", np.log(np.abs(cur_vars_df[new_feature]) + 1).values))
            for feature_name, feature in features_to_add:
                if feature_name in constants.all_cat_cols:
                    for i in range(feature.shape[1]):
                        cur_features_df[feature_name + f"_{i}"] = feature[:, i]
                else:
                    cur_features_df[feature_name] = feature
                all_r2 = []
                for seed in range(n_seeds):
                    # split on variants
                    train_variants, test_variants = train_test_split(cur_unique_vars, train_size=train_size, random_state=seed)
                    train_split, test_split = cur_features_df[cur_features_df.var_id.isin(train_variants)], cur_features_df[cur_features_df.var_id.isin(test_variants)]
                    Xs_train, Xs_test, ys_train, ys_test = train_split.iloc[:, 2:].values, test_split.iloc[:, 2:].values, train_split.predictand, test_split.predictand

                    # for categorical phenotypes use logit regression
                    model = LinearRegression() if phenotype_continuous else LogitRegression()
                    regr = model.fit(Xs_train, ys_train)
                    cur_r2 = r2_score(ys_test, regr.predict(Xs_test))
                    all_r2.append(cur_r2)
                cur_mean_r2 = np.mean(all_r2)
                # the feature has to improve mean r2 by at least half a percentage point
                if cur_mean_r2 > cur_best_mean_r2 + 0.005:
                    cur_best_mean_r2 = cur_mean_r2
                    cur_best_median_r2 = np.median(all_r2)
                    cur_best_std_r2 = np.std(all_r2)
                    cur_best_features.append(feature_name)
                else:
                    if feature_name in constants.all_cat_cols:
                        for i in range(feature.shape[1]):
                            cur_features_df.drop(columns=feature_name + f"_{i}", inplace=True)
                    else:
                        cur_features_df.drop(columns=feature_name, inplace=True)
        results.append((phenotype_gene, cur_best_features, cur_best_mean_r2, cur_best_median_r2, cur_best_std_r2, len(cur_vars_df), maf))
    return results

with Pool(100) as p:
    out = p.map(process_phenotype_gene, top_lr_results["phenotype/gene"])

results_df = pd.DataFrame(itertools.chain.from_iterable(out),
                          columns=["phenotype/gene", "features", "mean", "median", "std", "#variants", "maf"]
                         ).sort_values("mean", ascending=False).reset_index(drop=True)

results_df.to_pickle("regression_multi_feature.pkl")