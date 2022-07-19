import glob
import os

import pandas as pd

from pandas import isnull

from tqdm import tqdm
from multiprocessing import Pool

import constants
import myvariant

def fetch_variant(var_id):
    return mv.getvariant(var_id, assembly="hg38")


pp2_output_filename = constants.pp2_output_path

pp2_variants_df = pd.read_csv(pp2_output_filename, sep="\t")

pp2_variants_df.columns = list(map(lambda x: str.strip(x), pp2_variants_df.columns.tolist()))

pp2_variants_df.gene = pp2_variants_df["gene"].apply(str.strip)

all_adjusted_phenotypes = pd.read_pickle(constants.all_adjusted_phenotypes_path)

genes_for_phenotypes = set(all_adjusted_phenotypes.explode("genes").genes.unique())

# add mapping between genes according to aliases
gene_name_mapping = {"TLR9": 'AC097637.1', "GFUS": "TSTA3", "POLR1H": "ZNRD1"}
genes_for_phenotypes = genes_for_phenotypes.union(list(gene_name_mapping.keys()))

# take only the variants in genes of interest
pp2_variants_df = pp2_variants_df[pp2_variants_df.gene.isin(genes_for_phenotypes)]

pp2_variants_df = pp2_variants_df[constants.pp2_id_cols + constants.pp2_feature_cols_cat + constants.pp2_feature_cols_cont]

humdiv_pp2 = pd.read_csv(constants.pp2_humdiv_output_path, sep="\t")

humdiv_pp2.columns = list(map(lambda x: str.strip(x), humdiv_pp2.columns.tolist()))

humdiv_pp2.gene = humdiv_pp2["gene"].apply(str.strip)

pp2_output_id_cols = ["pos", "aa2", "gene"]

pp2_variants_df = pp2_variants_df.merge(humdiv_pp2[pp2_output_id_cols + constants.pp2_pred_cols], on=pp2_output_id_cols, how="left")

humvar_pp2 = pd.read_csv(constants.pp2_humvar_output_path, sep="\t")

humvar_pp2.columns = list(map(lambda x: str.strip(x), humvar_pp2.columns.tolist()))

humvar_pp2.gene = humvar_pp2["gene"].apply(str.strip)

pp2_variants_df = pp2_variants_df.merge(humvar_pp2[pp2_output_id_cols + constants.pp2_pred_cols], on=pp2_output_id_cols, how="left")

pp2_variants_df.pos = pp2_variants_df.pos.astype(int)
for col_name in ["aa1", "aa2", "acc", "chr_pos", "gene", "dbMAF"]:
    pp2_variants_df[col_name] = pp2_variants_df[col_name].apply(str.strip)
pp2_variants_df.dbMAF = pd.to_numeric(pp2_variants_df.dbMAF, errors="coerce")

for col_name in constants.pp2_feature_cols_cat:
    pp2_variants_df[col_name] = pp2_variants_df[col_name].astype(str).apply(str.strip).astype("category").cat.codes

# add another input to indicate that and set NAs to zero
for field_name in constants.pp2_feature_cols_cont:
    pp2_variants_df[field_name] = pd.DataFrame(pd.to_numeric(pp2_variants_df[field_name].astype(str).apply(str.strip)))
    pp2_variants_df[f"{field_name}_na"] = 0
    pp2_variants_df[f"{field_name}_na"][pp2_variants_df[field_name].isna()] = 1
    pp2_variants_df[field_name].fillna(0, inplace=True)

pp2_variants_df["var_id"] = pp2_variants_df.chr_pos.apply(str.strip) + pp2_variants_df.nt1.apply(str.strip) + "/" + pp2_variants_df.nt2.apply(str.strip)

pp2_variants_df.gene = pp2_variants_df.gene.apply(lambda x: gene_name_mapping[x] if x in gene_name_mapping else x)

pp2_variants_df = pp2_variants_df.drop(columns=["chr_pos", "nt1", "nt2"])
pp2_variants_df.drop_duplicates("var_id", inplace=True)

var_ids = pp2_variants_df.var_id

# get EVE scores
eve_files = glob.glob(f"{constants.eve_scores_dir}/*")

eve_genes = set(list(map(lambda x: x.rsplit("/", 1)[1].split("_", 1)[0], eve_files)))

# include EVE
eve_df = None
for gene in pp2_variants_df.gene.unique():
    gene_eve_path = f"{constants.eve_scores_dir}/{gene}_HUMAN.csv"
    if os.path.isfile(gene_eve_path):
        cur_eve_df = pd.read_csv(gene_eve_path)[["wt_aa", "position", "mt_aa", "EVE_scores_ASM", "uncertainty_ASM"]]
        cur_eve_df["gene"] = gene
        if eve_df is None:
            eve_df = cur_eve_df
        else:
            eve_df = pd.concat((eve_df, cur_eve_df), axis=0)

eve_df = eve_df[~(eve_df.EVE_scores_ASM.isna() & eve_df.uncertainty_ASM.isna())]

eve_df.rename(columns={"wt_aa": "aa1", "mt_aa": "aa2", "position": "pos"}, inplace=True)

pp2_variants_df = pp2_variants_df.merge(eve_df, on=["aa1", "aa2", "pos", "gene"], how="left")

# add cols to indicate NAs and set NAs to zero
for field_name in constants.eve_cols:
    pp2_variants_df[f"{field_name}_na"] = 0
    pp2_variants_df[f"{field_name}_na"][pp2_variants_df[field_name].isna()] = 1
    pp2_variants_df[field_name].fillna(0, inplace=True)

# add "_na" fields for the remaining continuous columns
for field_name in constants.pp2_pred_cols_both:
    pp2_variants_df[field_name] = pd.to_numeric(pp2_variants_df[field_name], errors="coerce")
    pp2_variants_df[f"{field_name}_na"] = 0
    pp2_variants_df[f"{field_name}_na"][pp2_variants_df[field_name].isna()] = 1
    pp2_variants_df[field_name].fillna(0, inplace=True)

pp2_variants_df.drop_duplicates("var_id", inplace=True)
# save all the variant features
pp2_variants_df.to_pickle(constants.variant_features_all)

# save vcf input
chromosome = pp2_variants_df.var_id.apply(lambda x: x[3:].split(":")[0])
pos = pp2_variants_df.var_id.apply(lambda x: x.split(":")[1][:-3])
ref = pp2_variants_df.var_id.apply(lambda x: x[-3])
alt = pp2_variants_df.var_id.apply(lambda x: x[-1])

vcf_input = pd.concat((chromosome, pos, pd.Series(["."]*len(pos)), ref, alt), axis=1)

vcf_input["QUAL"] = "."
vcf_input["FILTER"] = "PASS"
vcf_input["INFO"] = "."

vcf_input.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']

vcf_input["ID"] = "."

vcf_input.to_csv(constants.vcf_input_path, sep="\t", index=False)