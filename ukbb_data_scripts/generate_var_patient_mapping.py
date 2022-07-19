import pandas as pd
import json
import glob
from tqdm import tqdm

import constants

def create_chr_desc(x):
    chromosome, pos, ref, alt = x.split("-")
    return f"chr{chromosome}:{pos}{ref}/{alt}"

set_variants = set(pd.read_pickle(constants.variant_features_all).var_id.tolist())
mt_dir = constants.ukbb_var_patient_dir
mt_files = glob.glob(mt_dir + "/*")
variant_patient_mapping = {}
for mt_filename in tqdm(mt_files):
    mt_read = pd.read_csv(mt_filename)
    mt_read = mt_read[~mt_read.patient.isna()]
    mt_read.variant = mt_read.variant.apply(create_chr_desc).values
    mt_read = mt_read[mt_read.variant.isin(set_variants)]
    variants = mt_read.variant
    patients = mt_read.patient.apply(lambda x: list(map(int, str(x).split("|")))).values
    for variant, patient in zip(variants, patients):
        variant_patient_mapping[variant] = patient

# take at most 0.1 MAF
variant_patient_mapping_rare = dict()
for (key, value) in variant_patient_mapping.items():
    if len(value) < constants.maf_filter * constants.total_patients:
        variant_patient_mapping_rare[key] = value

with open(constants.var_patient_mapping_path, 'w') as fp:
    json.dump(variant_patient_mapping_rare, fp)
