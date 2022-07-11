import os
import glob

import pandas as pd
import numpy as np

from tqdm import tqdm
from sklearn.linear_model import LinearRegression

import constants
import helpers

# adjust for covariates
# adjust for sex, age, age^2, sex*age, sex*age^2, 10PCs
genebass_top_df = helpers.get_genebass_top()

pcs_df = pd.read_csv(constants.genetic_pcs_path)
age_df = pd.read_pickle(f"{constants.cleaned_adjust_phenotypes_dir}/{constants.age_phenocode}.pkl")
sex_df = pd.read_pickle(f"{constants.cleaned_adjust_phenotypes_dir}/{constants.sex_phenocode}.pkl")
statins_df = pd.read_pickle(f"{constants.cleaned_categorical_phenotypes_dir}/20003-statins.pkl").rename(
    columns={"value": "statins"})

sex_and_pcs_df = sex_df.merge(pcs_df, on="eid", how="left").drop(columns="assessment")

sex_and_pcs_df["pcs_given"] = ~sex_and_pcs_df.PC1.isna()

sex_and_pcs_df.rename(columns={"value": "sex"}, inplace=True)

adj_covariates = age_df.merge(sex_and_pcs_df, on="eid", how="left")

adj_covariates.rename(columns={"value": "age"}, inplace=True)

adj_covariates["age^2"] = adj_covariates.age * adj_covariates.age
adj_covariates["sex*age"] = adj_covariates.age * adj_covariates.sex
adj_covariates["sex*age^2"] = adj_covariates["age^2"] * adj_covariates.sex

# set pcs to the pcs of the person who has median phenotype value where not given
def fill_missing_pcs(phenotype_df):
    vals_for_median = phenotype_df[phenotype_df.pcs_given].value
    if len(vals_for_median) % 2 == 0:
        vals_for_median = vals_for_median[:-1]
    median_row = phenotype_df[np.isclose(phenotype_df.value, vals_for_median.median()) &
                             phenotype_df.pcs_given].iloc[0]
    for i in range(1, 11):
        phenotype_df.loc[phenotype_df.pcs_given == False, f"PC{i}"] = median_row[f"PC{i}"]

# adjust regular continuous variables with multiple measurements
phenotypes_files = glob.glob(f"{constants.cleaned_continuous_phenotypes_dir}/*")

adjustment_scores = []

for phenotype_path in tqdm(phenotypes_files):
    phenotype = phenotype_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    phenotype_df = pd.read_pickle(phenotype_path)
    phenotype_df = phenotype_df.merge(adj_covariates, on=["eid", "assessment"])
    if phenotype in genebass_top_df[genebass_top_df.Description.isin(constants.phenotypes_adjusted_for_statins)].Phenotype.unique():
        phenotype_df = phenotype_df.merge(statins_df, on="eid")
    if phenotype == "40007":
        phenotype_df = phenotype_df.drop(columns=["age", "age^2", "sex*age", "sex*age^2"])
    
    fill_missing_pcs(phenotype_df)
    reg = LinearRegression().fit(phenotype_df.drop(columns=["eid", "value", "assessment", "pcs_given"]), phenotype_df.value)
    r2 = reg.score(phenotype_df.drop(columns=["eid", "value", "assessment", "pcs_given"]), phenotype_df.value)
    save_path = f"{constants.adjusted_continuous_phenotypes_dir}/{phenotype}.pkl"
    adjustment_scores.append((phenotype, r2, "continuous", save_path, phenotype_path))
    phenotype_df["adjusted"] = phenotype_df.value - reg.predict(phenotype_df.drop(columns=["eid", "value", "assessment", "pcs_given"]))
    phenotype_df[["eid", "adjusted", "assessment"]].to_pickle(save_path)
    

adj_covariates = adj_covariates[adj_covariates.assessment == '0']

# adjust continuous variables with single measurements
phenotypes_files = glob.glob(f"{constants.cleaned_continuous_single_assessment_phenotypes_dir}/*")
for phenotype_path in tqdm(phenotypes_files):
    phenotype = phenotype_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    phenotype_df = pd.read_pickle(phenotype_path)
    phenotype_df = phenotype_df.merge(adj_covariates, on="eid")
    phenotype_df.drop(columns="assessment", inplace=True)
    
    fill_missing_pcs(phenotype_df)
    reg = LinearRegression().fit(phenotype_df.drop(columns=["eid", "value", "pcs_given"]), phenotype_df.value)
    r2 = reg.score(phenotype_df.drop(columns=["eid", "value", "pcs_given"]), phenotype_df.value)
    save_path = f"{constants.adjusted_continuous_phenotypes_dir}/{phenotype}.pkl"
    adjustment_scores.append((phenotype, r2, "continuous", save_path, phenotype_path))
    phenotype_df["adjusted"] = phenotype_df.value - reg.predict(phenotype_df.drop(columns=["eid", "value", "pcs_given"]))
    phenotype_df[["eid", "adjusted"]].to_pickle(save_path)

adjustment_scores_cat = []
# go through binary variables
phenotypes_files = glob.glob(f"{constants.cleaned_categorical_phenotypes_dir}/*")
for phenotype_path in tqdm(phenotypes_files):
    phenotype = phenotype_path.rsplit("/", 1)[-1].rsplit(".", 1)[0]
    phenotype_df = pd.read_pickle(phenotype_path)
    save_path = f"{constants.adjusted_categorical_phenotypes_dir}/{phenotype}.pkl"
    adjustment_scores_cat.append((phenotype, 0, "categorical", save_path, phenotype_path))
    phenotype_df["adjusted"] = phenotype_df.value
    phenotype_df[["eid", "adjusted"]].to_pickle(save_path)

# make phenotype-name and phenotype-gene mapping
phenocode_phenotype_map = dict(zip(genebass_top_df.Phenotype, genebass_top_df.Description))

phenocode_genes_map = dict(genebass_top_df.groupby("Phenotype").Gene.apply(list))

results_df = pd.DataFrame(adjustment_scores, columns=["phenocode", "adjustment_r2", "type", "adjusted_path", "raw_path"]).sort_values("adjustment_r2", ascending=False)

results_df["phenotype"] = results_df.phenocode.apply(lambda x: phenocode_phenotype_map[x])

results_df = results_df[results_df.phenotype != 'Year of birth']

results_df["genes"] = results_df.phenocode.apply(lambda x: phenocode_genes_map[x]) 



# categorical results processing

results_cat_df = pd.DataFrame(adjustment_scores_cat, columns=["phenocode", "adjustment_r2", "type", "adjusted_path", "raw_path"]).sort_values("adjustment_r2", ascending=False)

results_cat_df["phenotype"] = results_cat_df.phenocode
results_cat_df.loc[results_cat_df.phenocode == constants.wears_glasses_phenocode, "phenotype"] = \
    "Wears glasses or contact lenses"

results_cat_df["genes"] = [[] for i in range(len(results_cat_df))]

results_cat_df.loc[results_cat_df.phenocode.isin(phenocode_genes_map), "genes"] = \
    results_cat_df.loc[results_cat_df.phenocode.isin(phenocode_genes_map), "phenocode"].apply(lambda x: phenocode_genes_map[x])

# deal with the wears glassess trait by assigning the same genes as for the age when started wearing glassess
results_cat_df.loc[results_cat_df.phenocode == constants.wears_glasses_phenocode, "genes"] = \
    results_cat_df.loc[results_cat_df.phenocode == constants.wears_glasses_phenocode, "phenocode"].apply(lambda x: phenocode_genes_map['2217'])


# use the same genes as for LDL direct for the statins
results_cat_df.loc[results_cat_df.phenocode.apply(lambda x: x.split("-")[1] == "statins" if "-" in x else False), "genes"] = \
    results_cat_df[results_cat_df.phenocode.apply(lambda x: x.split("-")[1] == "statins" if "-" in x else False)].phenocode.apply(lambda x: phenocode_genes_map['30780'])


cat_only_genebass = genebass_top_df[(genebass_top_df["Trait type"] == "categorical")]
cat_only_genebass["joint_phenocode"] = cat_only_genebass.Phenotype + "-" + cat_only_genebass.analysis_id.apply(lambda x: x.split("-")[-2])
complex_phenocode_genes_map = dict(cat_only_genebass.groupby("joint_phenocode").Gene.apply(list))

results_cat_df.loc[results_cat_df.phenocode.apply(lambda x: "-" in x and x in complex_phenocode_genes_map), "genes"] = \
    results_cat_df.loc[results_cat_df.phenocode.apply(lambda x: "-" in x and x in complex_phenocode_genes_map), "phenocode"].apply(lambda x: complex_phenocode_genes_map[x])


pd.concat((results_df, results_cat_df), axis=0).reset_index(drop=True).to_pickle(constants.all_adjusted_phenotypes_path)
