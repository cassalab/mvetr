repo_dir = "/net/home/aszalata/repos/lipids"
basket_path = "/net/ukbb/datasets/49916/ukb49916.csv"
genebass_top_hist_path = f"{repo_dir}/data/gene-phewas-exomes_top-hits_2022_06_24_15_39_05.csv"
ukbb_extracts_output_dir = f"{repo_dir}/data/ukbb_extract"
genetic_pcs_path = f"{repo_dir}/data/PC_1_10_UKB.csv"
medication_coding_path = f"{repo_dir}/data/coding4.tsv"
eve_scores_dir = f"{repo_dir}/notebooks/gene_exploration/variant_files"
vcf_input_path = f"{repo_dir}/dbnsfp_scripts/my_input.vcf"
vcf_output_path = f"{repo_dir}/dbnsfp_scripts/my_output_all.vcf"

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

dbnsfp_categorical_cols = ["IMPACT", "FLAGS", "BayesDel_addAF_pred", "BayesDel_noAF_pred", "ClinPred_pred", "FATHMM_pred", "LRT_pred", "M-CAP_pred",
                          "MetaLR_pred", "MetaSVM_pred", "PrimateAI_pred", "fathmm-MKL_coding_group", "fathmm-MKL_coding_pred",
                          "fathmm-XF_coding_pred", "gnomAD_exomes_flag", "gnomAD_genomes_flag", "Chrom", "STRAND"]
dbnsfp_continuous_cols = ['1000Gp3_AF',
 'ALSPAC_AF',
 'BayesDel_addAF_rankscore',
 'BayesDel_addAF_score',
 'BayesDel_noAF_rankscore',
 'BayesDel_noAF_score',
 'CADD_phred',
 'CADD_phred_hg19',
 'CADD_raw',
 'CADD_raw_hg19',
 'CADD_raw_rankscore',
 'CADD_raw_rankscore_hg19',
 'ClinPred_rankscore',
 'ClinPred_score',
 'DANN_rankscore',
 'DANN_score',
 'DEOGEN2_rankscore',
 'DEOGEN2_score',
 'ESP6500_AA_AF',
 'ESP6500_EA_AF',
 'Eigen-PC-phred_coding',
 'Eigen-PC-raw_coding',
 'Eigen-PC-raw_coding_rankscore',
 'Eigen-phred_coding',
 'Eigen-raw_coding',
 'Eigen-raw_coding_rankscore',
 'ExAC_AF',
 'ExAC_nonTCGA_AF',
 'ExAC_nonpsych_AF',
 'FATHMM_converted_rankscore',
 'FATHMM_score',
 'GERP++_NR',
 'GERP++_RS',
 'GERP++_RS_rankscore',
 'GM12878_confidence_value',
 'GM12878_fitCons_rankscore',
 'GM12878_fitCons_score',
 'GenoCanyon_rankscore',
 'GenoCanyon_score',
 'H1-hESC_confidence_value',
 'H1-hESC_fitCons_rankscore',
 'H1-hESC_fitCons_score',
 'HUVEC_confidence_value',
 'HUVEC_fitCons_rankscore',
 'HUVEC_fitCons_score',
 'LINSIGHT',
 'LINSIGHT_rankscore',
 'LIST-S2_rankscore',
 'LIST-S2_score',
 'LRT_Omega',
 'LRT_converted_rankscore',
 'LRT_score',
 'M-CAP_rankscore',
 'M-CAP_score',
 'MPC_rankscore',
 'MPC_score',
 'MVP_rankscore',
 'MVP_score',
 'MetaLR_rankscore',
 'MetaLR_score',
 'MetaRNN_rankscore',
 'MetaRNN_score',
 'MetaSVM_rankscore',
 'MetaSVM_score',
 'MutPred_rankscore',
 'MutPred_score',
 'MutationAssessor_rankscore',
 'MutationAssessor_score',
 'MutationTaster_converted_rankscore',
 'MutationTaster_score',
 'PROVEAN_converted_rankscore',
 'PROVEAN_score',
 'Polyphen2_HDIV_rankscore',
 'Polyphen2_HDIV_score',
 'Polyphen2_HVAR_rankscore',
 'Polyphen2_HVAR_score',
 'PrimateAI_rankscore',
 'PrimateAI_score',
 'REVEL_rankscore',
 'REVEL_score',
 'Reliability_index',
 'SIFT4G_converted_rankscore',
 'SIFT4G_score',
 'SIFT_converted_rankscore',
 'SIFT_score',
 'SiPhy_29way_logOdds',
 'SiPhy_29way_logOdds_rankscore',
 'TSL',
 'TWINSUK_AC',
 'TWINSUK_AF',
 'UK10K_AC',
 'UK10K_AF',
 'VEST4_rankscore',
 'VEST4_score',
 'aapos',
 'bStatistic',
 'bStatistic_converted_rankscore',
 'codon_degeneracy',
 'codonpos',
 'fathmm-MKL_coding_rankscore',
 'fathmm-MKL_coding_score',
 'fathmm-XF_coding_rankscore',
 'fathmm-XF_coding_score',
 'gnomAD_exomes_AF',
 'gnomAD_exomes_POPMAX_AF',
 'gnomAD_exomes_POPMAX_nhomalt',
 'gnomAD_exomes_controls_AF',
 'gnomAD_exomes_controls_POPMAX_AF',
 'gnomAD_exomes_controls_POPMAX_nhomalt',
 'gnomAD_exomes_controls_nhomalt',
 'gnomAD_exomes_nhomalt',
 'gnomAD_exomes_non_cancer_AF',
 'gnomAD_exomes_non_cancer_POPMAX_AF',
 'gnomAD_exomes_non_cancer_POPMAX_nhomalt',
 'gnomAD_exomes_non_cancer_nhomalt',
 'gnomAD_exomes_non_neuro_AF',
 'gnomAD_exomes_non_neuro_POPMAX_AF',
 'gnomAD_exomes_non_neuro_POPMAX_nhomalt',
 'gnomAD_exomes_non_neuro_nhomalt',
 'gnomAD_exomes_non_topmed_AF',
 'gnomAD_exomes_non_topmed_nhomalt',
 'gnomAD_genomes_AF',
 'gnomAD_genomes_POPMAX_AF',
 'gnomAD_genomes_POPMAX_nhomalt',
 'gnomAD_genomes_controls_and_biobanks_AF',
 'gnomAD_genomes_controls_and_biobanks_nhomalt',
 'gnomAD_genomes_nhomalt',
 'gnomAD_genomes_non_cancer_AF',
 'gnomAD_genomes_non_cancer_nhomalt',
 'gnomAD_genomes_non_neuro_AF',
 'gnomAD_genomes_non_topmed_nhomalt',
 'integrated_confidence_value',
 'integrated_fitCons_rankscore',
 'integrated_fitCons_score',
 'phastCons100way_vertebrate',
 'phastCons100way_vertebrate_rankscore',
 'phastCons17way_primate',
 'phastCons17way_primate_rankscore',
 'phastCons30way_mammalian',
 'phastCons30way_mammalian_rankscore',
 'phyloP100way_vertebrate',
 'phyloP100way_vertebrate_rankscore',
 'phyloP17way_primate',
 'phyloP17way_primate_rankscore',
 'phyloP30way_mammalian',
 'phyloP30way_mammalian_rankscore',
 'canonical_found']

varity_cols = ['VARITY_ER', 'VARITY_ER_LOO', 'provean_score', 'sift_score',
       'evm_epistatic_score', 'GERP_RS', 'blosum100', 'in_domain', 'asa_mean',
       'aa_psipred_E', 'aa_psipred_H', 'aa_psipred_C', 'bsa_max', 'h_bond_max',
       'salt_bridge_max', 'disulfide_bond_max', 'covelent_bond_max',
       'solv_ne_abs_max', 'mw_delta', 'pka_delta', 'pkb_delta', 'pi_delta',
       'hi_delta', 'pbr_delta', 'avbr_delta', 'vadw_delta', 'asa_delta',
       'cyclic_delta', 'charge_delta', 'positive_delta', 'negative_delta',
       'hydrophobic_delta', 'polar_delta', 'ionizable_delta', 'aromatic_delta',
       'aliphatic_delta', 'hbond_delta', 'sulfur_delta', 'essential_delta',
       'size_delta']

all_cat_cols = pp2_feature_cols_cat + dbnsfp_categorical_cols
all_cont_cols = pp2_feature_cols_cont + eve_cols + pp2_pred_cols_both + dbnsfp_continuous_cols + varity_cols
