repo_dir = "/net/home/aszalata/repos/lipids"
basket_path = "/net/ukbb/datasets/49916/ukb49916.csv"
genebass_top_hist_path = f"{repo_dir}/data/gene-phewas-exomes_top-hits_2022_06_24_15_39_05.csv"
ukbb_extracts_output_dir = f"{repo_dir}/data/ukbb_extract"
genetic_pcs_path = f"{repo_dir}/data/PC_1_10_UKB.csv"
medication_coding_path = f"{repo_dir}/data/coding4.tsv"
eve_scores_dir = f"{repo_dir}/notebooks/gene_exploration/variant_files"
varity_scores_path = f"{repo_dir}/notebooks/data/datasets/varity_all_predictions.txt"
vcf_input_path = f"{repo_dir}/dbnsfp_scripts/my_input.vcf"
vcf_input_path = f"{repo_dir}/dbnsfp_scripts/my_output_all.vcf"

vcfs_dir = "/net/ukbb/datasets/45878/hail_extracts/non_annotated_vcfs"
pp2_output_path = "/net/data/artur/pph2/afold2/mysnps.pph.output"
pp2_humdiv_output_path = "/net/data/artur/pph2/afold2/mysnps.humdivmz.output"
pp2_humvar_output_path = "/net/data/artur/pph2/afold2/mysnps.humvarmz.output"
# all of the variant features we use
variant_features_all = f"{repo_dir}/data/variant_features_all.pkl"

# path to the df with the missense variants ids from pp2
pp2_variants_path = f"{repo_dir}/data/ukbb_extract/pp2_variants.pkl"

cleaned_adjust_phenotypes_dir = f"{ukbb_extracts_output_dir}/cleaned_phenotypes/adjust"
cleaned_continuous_phenotypes_dir = f"{ukbb_extracts_output_dir}/cleaned_phenotypes/continuous"
cleaned_continuous_single_assessment_phenotypes_dir = f"{ukbb_extracts_output_dir}/cleaned_phenotypes/continuous_single_assessment"
cleaned_categorical_phenotypes_dir = f"{ukbb_extracts_output_dir}/cleaned_phenotypes/categorical_binary"

all_adjusted_phenotypes_path = f"{ukbb_extracts_output_dir}/adjusted_phenotypes/all_adjusted_df.pkl"

adjusted_continuous_phenotypes_dir = f"{ukbb_extracts_output_dir}/adjusted_phenotypes/continuous"
adjusted_categorical_phenotypes_dir = f"{ukbb_extracts_output_dir}/adjusted_phenotypes/categorical"

age_phenocode = "21003"
sex_phenocode = "31"
medication_phenocode = "20003"
adjustment_phenocodes = [age_phenocode, sex_phenocode]

statins_names = ["atorvastatin", "fluvastatin", "lovastatin", "pitavastatin", "pravastatin", "rosuvastatin", "simvastatin", "lipitor", "altoprev", "lescol", "livalo", "zypitamag", "pravachol", "crestor", "ezallor", "zocor"]


# subcategories of genebass phenotypes
continuous_phenotypes_other = [
'Exposure to tobacco smoke at home', # take only 0 and greater, continuous
'Age when periods started (menarche)', # only where greater than 0, continuous
'Age at menopause (last menstrual period)', # only answers greater than 0, continuous
'Age started wearing glasses or contact lenses', # continuous only for those who wear reject nan and smaller than 0
'Recent changes in speed/amount of moving or speaking', # only answers containing 1-4, categorical
'Recent lack of interest or pleasure in doing things', # only 1-4, categorical
'Frequency of failure to fulfil normal expectations due to drinking alcohol in last year', # only 1-5, categorical
'Recent trouble concentrating on things', # 1-4 categorical
'Recent feelings of depression', # 1-4 cat
'Recent feelings of foreboding', # greater than 0, cat
'Frequency of memory loss due to drinking alcohol in last year', # cat, only greater than 0
'Smoking/smokers in household', # cat, 0 or greater
'Ever had known person concerned about, or recommend reduction of, alcohol consumption', # cat, only greater than 0
'Recent easy annoyance or irritability', # cat, only greater than 0
'Frequency of feeling guilt or remorse after drinking alcohol in last year', # cat, only greater than 0
'Recent poor appetite or overeating', # cat only greater than 0
'Recent inability to stop or control worrying'] # cat only greater than 0

categorical_multiclass_phenocodes = [
    # in categorical, take categories according to top phenotype-gene, make binary
# categorical:
    "6152",
    "6177",
    "6150",
    "6148",
    "6154",
]

categorical_text_multiclass_phenocodes = [
# look into genebass description to find particular associations
# in categorical - collect all of the patient-procedure pairs and use categories instead of names? binary then
# create additional category - any of the disease
# assume not present means negative class
# binary:
"41200",
"20004",
"20002",
"41210",
"20001", 
"20003",
"20544"
]

# special binary for which we use the genebass association corresponding to 'Age started wearing glasses or contact lenses'
wears_glasses_phenocode = "2207"

categorical_binary_phenocodes = [
    "2443", # binary, take only 0 and 1
    "2463", # binary, take only 0 and 1
]

phenotypes_adjusted_for_statins = [
    "Cholesterol", "LDL direct", "Apolipoprotein B", "HDL cholesterol", "Apolipoprotein A"
]

# PP2 constants
pp2_id_cols = ["pos", "aa1", "aa2", "acc", "chr_pos", "gene", 'nt1', 'nt2', "dbMAF"]
pp2_feature_cols_cat = ["dgn", "prediction", "based_on", "effect", "site", "region", "JXc", "CpG", "str", "dref", "SecStr", "MapReg"]
pp2_feature_cols_cont = ["PHAT", "dScore", "JXmin", "trv", "phylop", "Score1", "Score2", "Nobs", "Nseqs", "Nsubs",
                     "Nvars", "Nres", "IdPmax", "IdPSNP", "IdQmax", "DistPmin", "DistPSNP", "DistQmin", "BaRE",
                     "Nstr", "Nfilt", "PDB_len", "PDB_pos", "PDB_idn", "dVol", "dProp", "NormASA", "B-fact", "H-bonds", "AveNSit", "MinDSit"]
pp2_pred_cols = ['pph2_prob', 'pph2_FPR', 'pph2_TPR', 'pph2_FDR']
# columns after merging humvar and humdiv pp2 predictions
pp2_pred_cols_both = ['pph2_prob_x', 'pph2_FPR_x', 'pph2_TPR_x', 'pph2_FDR_x', 'pph2_prob_y',
       'pph2_FPR_y', 'pph2_TPR_y', 'pph2_FDR_y']
eve_cols = ["EVE_scores_ASM", "uncertainty_ASM"]
varity_cols = ["VARITY_R", "VARITY_ER", "VARITY_R_LOO", "VARITY_ER_LOO"]
# excluding dbnsfp
all_cont_cols = pp2_feature_cols_cont + eve_cols + varity_cols + pp2_pred_cols_both

# dbnsfp columns without the "na" indicator columns
# mind that dbnsfp_raw_cols and pp2_feature_cols_cont have corresponding indicator cols with "_na" appendix
# pp2_feature_cols_cat are categorical columns for which we have to apply one hot encoding before using them
dbnsfp_raw_cols_path = f"{repo_dir}/data/ukbb_extract/dbnsfp_raw_cols.pkl"
