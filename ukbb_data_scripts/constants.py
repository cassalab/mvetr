repo_dir = "/net/home/aszalata/repos/lipids"
basket_path = "/net/ukbb/datasets/49916/ukb49916.csv"
genebass_top_hist_path = f"{repo_dir}/data/gene-phewas-exomes_top-hits_2022_06_24_15_39_05.csv"
ukbb_extracts_output_dir = f"{repo_dir}/data/ukbb_extract"
genetic_pcs_path = f"{repo_dir}/data/PC_1_10_UKB.csv"

age_phenocode = "21003"
sex_phenocode = "31"
adjustment_phenocodes = [age_phenocode, sex_phenocode]

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