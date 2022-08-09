import pandas as pd
import numpy as np

from ukbb_data_scripts import constants
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score
import global_constants

import itertools

from multiprocessing import Pool

import data_helpers

all_vars_phenotypes_df = pd.read_pickle(global_constants.top_regression_vars_phenotypes_path)

top_lr_results = pd.read_pickle("regression_results_w_multi_features_processed.pkl")


def process_phenotype_gene_pairs(args):
    train_size = 0.8
    n_seeds = 100
    phenotype_1 = args["phenotype_1"]
    phenotype_2 = args["phenotype_2"]
    gene_1 = args["gene_1"]
    gene_2 = args["gene_2"]
    maf_1 = args["maf_1"]
    maf_2 = args["maf_2"]
    features_1 = args["features_1"]
    features_2 = args["features_2"]
    phenotype_gene_1 = phenotype_1 + "/" + gene_1
    phenotype_gene_2 = phenotype_2 + "/" + gene_2
    
    results = []
    vars_df_1 = all_vars_phenotypes_df[(all_vars_phenotypes_df.phenotype == phenotype_1) & (all_vars_phenotypes_df.gene == gene_1) & (all_vars_phenotypes_df.highest_AF < maf_1)]
    vars_df_2 = all_vars_phenotypes_df[(all_vars_phenotypes_df.phenotype == phenotype_2) & (all_vars_phenotypes_df.gene == gene_2) & (all_vars_phenotypes_df.highest_AF < maf_2)]
    vars_df = pd.concat((vars_df_1, vars_df_2), axis=0)
    vars_df["phenotype/gene"] = vars_df["phenotype"] + "/" + vars_df["gene"]
    del vars_df_1
    del vars_df_2
    # try either features_1, features_2 or the union
    unique_vars = vars_df.var_id.unique()
    for features in [features_1, features_2, list(set(features_1).union(features_2))]:
        cur_vars_df = data_helpers.get_features_w_id_df(vars_df, features)
        for scale_outputs in [True, False]:
            all_r2_1 = []
            all_r2_2 = []
            for seed in range(n_seeds):
                # split on variants
                train_variants, test_variants = train_test_split(unique_vars, train_size=train_size, random_state=seed)
                train_split, test_split = cur_vars_df[cur_vars_df.var_id.isin(train_variants)], cur_vars_df[cur_vars_df.var_id.isin(test_variants)]
                Xs_train, Xs_test, ys_train, ys_test = train_split.iloc[:, 3:].values, test_split.iloc[:, 3:].values, train_split.predictand.values.reshape((-1, 1)), test_split.predictand.values.reshape((-1, 1))
                split_cond_1 = test_split["phenotype/gene"] == phenotype_gene_1
                split_cond_2 = test_split["phenotype/gene"] == phenotype_gene_2
                split_cond_1_train = train_split["phenotype/gene"] == phenotype_gene_1
                split_cond_2_train = train_split["phenotype/gene"] == phenotype_gene_2
                if split_cond_1_train.sum() < 2 or split_cond_2_train.sum() < 2 \
                        or split_cond_1.sum() < 2 or split_cond_1.sum() < 2:
                    continue
                
                if scale_outputs == True:
                    # separate scaling for each phenotype/gene
                    scaler_1 = StandardScaler().fit(ys_train[split_cond_1_train])
                    scaler_2 = StandardScaler().fit(ys_train[split_cond_2_train])
                    ys_train[split_cond_1_train] = scaler_1.transform(ys_train[split_cond_1_train])
                    ys_train[split_cond_2_train] = scaler_2.transform(ys_train[split_cond_2_train])
                    ys_test[split_cond_1] = scaler_1.transform(ys_test[split_cond_1])
                    ys_test[split_cond_2] = scaler_2.transform(ys_test[split_cond_2])
                else:
                    pass

                model = LinearRegression()
                regr = model.fit(Xs_train, ys_train)
                cur_r2_1 = r2_score(ys_test[split_cond_1], regr.predict(Xs_test[split_cond_1]))
                all_r2_1.append(cur_r2_1)
                cur_r2_2 = r2_score(ys_test[split_cond_2], regr.predict(Xs_test[split_cond_2]))
                all_r2_2.append(cur_r2_2)
            cur_mean_r2_1 = np.mean(all_r2_1)
            cur_median_r2_1 = np.median(all_r2_1)
            cur_std_r2_1 = np.std(all_r2_1)
            cur_mean_r2_2 = np.mean(all_r2_2)
            cur_median_r2_2 = np.median(all_r2_2)
            cur_std_r2_2 = np.std(all_r2_2)
            results.append((phenotype_gene_1, phenotype_gene_2, features, 
                            cur_mean_r2_1, cur_median_r2_1, cur_std_r2_1,
                            cur_mean_r2_2, cur_median_r2_2, cur_std_r2_2,
                            len(cur_vars_df), scale_outputs))
    return results


all_pairs = list(itertools.combinations(top_lr_results["phenotype/gene"], 2))
top_lr_results = top_lr_results[["phenotype", "gene", "maf", "features", "phenotype/gene"]].set_index("phenotype/gene")

tasks_list = []
for first_entry, second_entry in all_pairs:
    task = top_lr_results.loc[first_entry].add_suffix("_1").to_dict()
    task.update(top_lr_results.loc[second_entry].add_suffix("_2").to_dict())
    tasks_list.append(task)

with Pool(100) as p:
    out = p.map(process_phenotype_gene_pairs, tasks_list)

results_df = pd.DataFrame(itertools.chain.from_iterable(out),
                          columns=["phenotype/gene_1", "phenotype/gene_2",
                                   "features", "mean_1", "median_1", "std_1",
                                   "mean_2", "median_2", "std_2", "#variants", "z_outputs"]
                         ).sort_values("mean_1", ascending=False).reset_index(drop=True)

results_df.to_pickle("regression_phenotype_gene_pairs_results.pkl")
