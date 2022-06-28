import os
import glob

import pandas as pd

from tqdm import tqdm

import constants
import helpers


genebass_top_df = helpers.get_genebass_top()

# save pickled continuous phenotypes including the assessment
def clean_continuous_phenocodes(phenocode_list, output_dir):
    eid_df = pd.read_csv(os.path.join(constants.ukbb_extracts_output_dir, "eid.csv"))

    for phenocode in tqdm(phenocode_list):
        patient_phenotype_df = None
        phenotypes_files = glob.glob(constants.ukbb_extracts_output_dir + f"/{phenocode}-*.csv")
        for filename in phenotypes_files:
            version = filename.rsplit("/", 1)[1].split("-")[1].split(".")[0]
            df = pd.read_csv(filename, skip_blank_lines=False)
            df = pd.concat((eid_df, df), axis=1)
            df.dropna(inplace=True)
            df["assessment"] = version
            df.rename(columns={df.columns[1]: "value"}, inplace=True)
            if patient_phenotype_df is None:
                patient_phenotype_df = df
            else:
                patient_phenotype_df = pd.concat((patient_phenotype_df, df))
        patient_phenotype_df.assessment = patient_phenotype_df.assessment.astype("category")
        patient_phenotype_df.to_pickle(f"{output_dir}/{phenocode}.pkl")


# save continuous phenotypes taking only nonnegative values and only first assessment
def clean_continuous_single_measurement_phenocodes(phenocode_list, output_dir):
    eid_df = pd.read_csv(os.path.join(constants.ukbb_extracts_output_dir, "eid.csv"))

    for phenocode in tqdm(phenocode_list):
        patient_phenotype_df = None
        phenotypes_files = glob.glob(constants.ukbb_extracts_output_dir + f"/{phenocode}-0.0.csv")
        for filename in phenotypes_files:
            df = pd.read_csv(filename, skip_blank_lines=False)
            df = pd.concat((eid_df, df), axis=1)
            df.dropna(inplace=True)
            df.rename(columns={df.columns[1]: "value"}, inplace=True)
            df = df[df.value >= 0]
            if patient_phenotype_df is None:
                patient_phenotype_df = df
            else:
                patient_phenotype_df = pd.concat((patient_phenotype_df, df))
        patient_phenotype_df.to_pickle(f"{output_dir}/{phenocode}.pkl")


# save binary categorical phenotype taking only nonnegative values
def clean_categorical_binary_phenocodes(phenocode_list, output_dir):
    eid_df = pd.read_csv(os.path.join(constants.ukbb_extracts_output_dir, "eid.csv"))

    for phenocode in tqdm(phenocode_list):
        patient_phenotype_df = None
        phenotypes_files = glob.glob(constants.ukbb_extracts_output_dir + f"/{phenocode}-0.0.csv")
        for filename in phenotypes_files:
            df = pd.read_csv(filename, skip_blank_lines=False)
            df = pd.concat((eid_df, df), axis=1)
            df.dropna(inplace=True)
            df.rename(columns={df.columns[1]: "value"}, inplace=True)
            df = df[df.value >= 0]
            if patient_phenotype_df is None:
                patient_phenotype_df = df
            else:
                patient_phenotype_df = pd.concat((patient_phenotype_df, df))
        patient_phenotype_df.to_pickle(f"{output_dir}/{phenocode}.pkl")
        

# save multiclass categorical phenotype splitting into binary per-class files, 
# taking only values different than -3 (meaning questionnnaire not answered)
def clean_categorical_multiclass_phenocodes(phenocode_list, output_dir):
    eid_df = pd.read_csv(os.path.join(constants.ukbb_extracts_output_dir, "eid.csv"))

    for phenocode in tqdm(phenocode_list):
        associated_codes = genebass_top_df[genebass_top_df.Phenotype ==  phenocode].analysis_id.apply(
            lambda x: x.split("-")[-2]).astype(int).unique()
        
        patient_phenotype_df = None
        phenotypes_files = glob.glob(constants.ukbb_extracts_output_dir + f"/{phenocode}-0.*.csv")
        for filename in phenotypes_files:
            df = pd.read_csv(filename, skip_blank_lines=False)
            df = pd.concat((eid_df, df), axis=1)
            df.dropna(inplace=True)
            df.rename(columns={df.columns[1]: "value"}, inplace=True)
            df = df[df.value != -3]
            df.loc[~df.value.isin(associated_codes), "value"] = -7
            if patient_phenotype_df is None:
                patient_phenotype_df = df
            else:
                patient_phenotype_df = pd.concat((patient_phenotype_df, df))
        
        for code in associated_codes:
            if code not in patient_phenotype_df.value.unique():
                continue
            cur_patient_phenotype_df = patient_phenotype_df.copy()
            cur_patient_phenotype_df.loc[cur_patient_phenotype_df.value != code, "value"] = 0
            cur_patient_phenotype_df.loc[cur_patient_phenotype_df.value == code, "value"] = 1
            cur_patient_phenotype_df = cur_patient_phenotype_df.sort_values("value", ascending=False).drop_duplicates(subset="eid")
            cur_patient_phenotype_df.to_pickle(f"{output_dir}/{phenocode}-{code}.pkl")


# save multiclass textual categorical phenotype splitting into binary per-class files, 
# assuming that all eids that don't appear belong to the negative class
def clean_categorical_text_multiclass_phenocodes(phenocode_list, output_dir):
    eid_df = pd.read_csv(os.path.join(constants.ukbb_extracts_output_dir, "eid.csv"))

    for phenocode in tqdm(phenocode_list):
        associated_codes = genebass_top_df[genebass_top_df.Phenotype ==  phenocode].analysis_id.apply(
            lambda x: x.split("-")[-2]).astype(str).unique()
        
        patient_phenotype_df = None
        phenotypes_files = glob.glob(constants.ukbb_extracts_output_dir + f"/{phenocode}-0.*.csv")
        for filename in phenotypes_files:
            df = pd.read_csv(filename, skip_blank_lines=False)
            df = pd.concat((eid_df, df), axis=1)
            df.rename(columns={df.columns[1]: "value"}, inplace=True)
            df.value.fillna(0)
            df.value = df.value.astype(str)
            df.loc[~df.value.isin(associated_codes), "value"] = '0'
            if patient_phenotype_df is None:
                patient_phenotype_df = df
            else:
                patient_phenotype_df = pd.concat((patient_phenotype_df, df))
        
        for code in associated_codes:
            if code not in patient_phenotype_df.value.unique():
                continue
            cur_patient_phenotype_df = patient_phenotype_df.copy()
            cur_patient_phenotype_df.loc[cur_patient_phenotype_df.value != code, "value"] = '0'
            cur_patient_phenotype_df.loc[cur_patient_phenotype_df.value == code, "value"] = '1'
            cur_patient_phenotype_df.value = cur_patient_phenotype_df.value.astype(float)
            cur_patient_phenotype_df = cur_patient_phenotype_df.sort_values("value", ascending=False).drop_duplicates(subset="eid")
            cur_patient_phenotype_df.to_pickle(f"{output_dir}/{phenocode}-{code}.pkl")


regular_continuous_phenocodes = genebass_top_df[(genebass_top_df["Trait type"] == "continuous") & 
                                                ~genebass_top_df.Description.isin(constants.continuous_phenotypes_other)
                                               ].Phenotype.unique()

clean_continuous_phenocodes(regular_continuous_phenocodes, f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/continuous")

clean_continuous_phenocodes(constants.adjustment_phenocodes, f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/adjust")

single_assessment_continuous_phenocodes = genebass_top_df[(genebass_top_df["Trait type"] == "continuous") & 
                                                genebass_top_df.Description.isin(constants.continuous_phenotypes_other)
                                               ].Phenotype.unique()

clean_continuous_single_measurement_phenocodes(single_assessment_continuous_phenocodes, 
                                               f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/continuous_single_assessment")

clean_categorical_binary_phenocodes(constants.categorical_binary_phenocodes + [constants.wears_glasses_phenocode], 
                                    f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/categorical_binary")

clean_categorical_multiclass_phenocodes(constants.categorical_multiclass_phenocodes, 
                                        f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/categorical_binary")

clean_categorical_text_multiclass_phenocodes(constants.categorical_text_multiclass_phenocodes, 
                                             f"{constants.ukbb_extracts_output_dir}/cleaned_phenotypes/categorical_binary")