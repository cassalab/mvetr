import pandas as pd
import numpy as np
import itertools

from tqdm import tqdm
from multiprocessing import Pool
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

from ukbb_data_scripts import constants

from model_helpers import LogitRegression


all_adjusted_phenotypes = pd.read_pickle(constants.all_adjusted_phenotypes_path)

variant_features_all = pd.read_pickle(constants.variant_features_all)

def process_phenotype(phenotype):
    results = []
    phenotype_continuous = all_adjusted_phenotypes[all_adjusted_phenotypes.phenotype == phenotype].type.iloc[0] == "continuous"
    if phenotype_continuous:
        predictand = "median"
    else:
        predictand = "mean"
    vars_df = pd.read_pickle(all_adjusted_phenotypes[all_adjusted_phenotypes.phenotype == phenotype].var_phenotype_path.iloc[0])
    for gene in vars_df.gene.unique():
        cur_vars_df_no_maf = vars_df[vars_df.gene == gene][[predictand, "var_id"]].merge(
            variant_features_all[["var_id"] + constants.all_cat_cols + constants.all_cont_cols + list(map(lambda x: x + "_na", constants.all_cont_cols))], how="left")
        for maf in [0.1, 0.01, 0.001, 0.0001]:
            cur_vars_df = cur_vars_df_no_maf[cur_vars_df_no_maf.highest_AF < maf]
            if len(cur_vars_df) < 50:
                continue
            for extended_features in [True, False]:
                for feature in constants.all_cat_cols + constants.all_cont_cols:
                    train_size = 0.8
                    n_seeds = 100
                    best_r2 = -99999
                    all_r2 = []
                    r2_per_gene = []
                    for seed in range(n_seeds):
                        train_indices, test_indices = train_test_split(range(len(cur_vars_df)), train_size=train_size, random_state=seed)

                        if feature in constants.all_cat_cols:
                            pred_features = pd.get_dummies(cur_vars_df[feature]).values
                        else:
                            base_columns = cur_vars_df[[feature] + [feature + "_na"]]
                            if extended_features:
                                pred_features = pd.concat((base_columns, 
                                                           cur_vars_df[feature]**2,
                                                          np.log(np.abs(cur_vars_df[feature]) + 1)), axis=1).values
                            else:
                                pred_features = base_columns.values
                        Xs, ys = pred_features, cur_vars_df[predictand].values
                        Xs_train, Xs_test, ys_train, ys_test = Xs[train_indices], Xs[test_indices], ys[train_indices], ys[test_indices]

                        # for categorical phenotypes use logit regression
                        model = LinearRegression() if phenotype_continuous else LogitRegression()
                        regr = model.fit(Xs_train, ys_train)
                        cur_r2 = r2_score(ys_test, regr.predict(Xs_test))
                        all_r2.append(cur_r2)
                        if cur_r2 > best_r2:
                            best_r2 = cur_r2
                            best_regr = regr
                    results.append((phenotype + "/" + gene, feature, np.mean(all_r2), np.median(all_r2), np.std(all_r2), len(cur_vars_df), maf, extended_features))
    return results

with Pool(100) as p:
    out = p.map(process_phenotype, all_adjusted_phenotypes.phenotype.unique())

results_df = pd.DataFrame(itertools.chain.from_iterable(out),
                          columns=["phenotype/gene", "feature", "mean", "median", "std", "#variants", "maf", "extended_features"]
                         ).sort_values("median", ascending=False).reset_index(drop=True)

results_df.to_pickle("regression_results_extended.pkl")