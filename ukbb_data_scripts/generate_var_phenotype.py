import pandas as pd
import numpy as np
import json
import glob
from tqdm import tqdm

import constants


with open(constants.var_patient_mapping_path, 'r') as fp:
    variant_patient_mapping_rare = json.load(fp)
    
variant_freq = dict()
for (key, value) in variant_patient_mapping_rare.items():
    variant_freq[key] = len(value) / constants.total_patients

variant_features_all = pd.read_pickle(constants.variant_features_all)

variant_features_all = variant_features_all[variant_features_all.var_id.isin(variant_freq)]

variant_features_all["ukbb_AF"] = variant_features_all.var_id.apply(lambda x: variant_freq[x])

af_columns = list(filter(lambda x: x[-3:] == "_AF", variant_features_all.columns.tolist()))

# remove variants with MAF higher than 0.1
for af_col in af_columns:
    variant_features_all = variant_features_all[variant_features_all[af_col].apply(lambda x: x < 0.1)]

# create a column that gives the highest AF observed
variant_features_all["highest_AF"] = variant_features_all[af_columns].max(axis=1)
variant_features_all["highest_AF_na"] = 0

variant_features_all.to_pickle(constants.variant_features_all)

gene_var_mapping = variant_features_all.groupby("gene").var_id.apply(list).to_dict()

all_adjusted_phenotypes = pd.read_pickle(constants.all_adjusted_phenotypes_path)

# create dataframes per phenotype with variant-gene-median-mean-std(if more than 1 else None)-median/phenotype_std-median/std_of_medians
save_paths = []
for row in tqdm(all_adjusted_phenotypes.iterrows(), total=len(all_adjusted_phenotypes)):
    new_df = []
    phenotype = row[1].phenotype
    adjusted = pd.read_pickle(row[1].adjusted_path)
    genes = row[1].genes
    
    phenotype_std = adjusted.adjusted.std()
    for gene in genes:
        gene_vars_values = []
        if gene not in gene_var_mapping:
            continue
        gene_vars = gene_var_mapping[gene]
        for var in gene_vars:
            patient_ids = variant_patient_mapping_rare[var]
            values = adjusted[adjusted.eid.isin(patient_ids)].adjusted.values
            if len(values) == 0:
                continue
            gene_vars_values.append((var, values))
        if len(gene_vars_values) == 0:
            continue
        gene_std = np.std(list(map(lambda x: np.median(x[1]), gene_vars_values)))
        for var, values in gene_vars_values:
            new_df.append((var, gene, np.median(values), np.mean(values),
                           np.std(values) if len(values) > 1 else None,
                          np.median(values)/phenotype_std, np.median(values)/gene_std))
    new_df = pd.DataFrame(new_df, columns=["var_id", "gene", "median", "mean", "std",
                              "median/phenotype_std", "median/gene_medians_std"])
    save_path = f"{constants.variant_phenotype_output_dir}/{phenotype.replace("/", "_")}.pkl"
    save_paths.append(save_path)
    new_df.to_pickle(save_path)

all_adjusted_phenotypes["var_phenotype_path"] = save_paths

all_adjusted_phenotypes.to_pickle(constants.all_adjusted_phenotypes_path)