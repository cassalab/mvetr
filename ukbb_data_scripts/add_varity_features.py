import glob
import re

import pandas as pd
from tqdm import tqdm

from selenium import webdriver

import constants


regex = r"\{(.*?)\}"
DRIVER_PATH = '/usr/local/bin/chromedriver'
options = webdriver.ChromeOptions()
options.add_argument("headless")
driver = webdriver.Chrome(executable_path=DRIVER_PATH, options=options)

all_features_df = pd.read_pickle(constants.variant_features_all)
for uniprot_id in tqdm(all_features_df.acc.unique().values):
    driver.get(f'http://varity.varianteffect.org/wsgi/varity_wsgi.wsgi?queryflag=9;sessionid={uniprot_id};search_type=uniprot_id;score_name=VARITY_R;from_pos=1;to_pos=99999999999;callback=callback1')

    matches = re.finditer(regex, driver.page_source, re.MULTILINE | re.DOTALL)

    for matchNum, match in enumerate(matches):
        for groupNum in range(0, len(match.groups())):
            download_path = eval("{" + match.group(1) + "}")["content"]

    driver.get(f"http://varity.varianteffect.org/output/{download_path}.csv")

varity_paths = glob.glob("*.csv")

cur_aa = pd.DataFrame(all_features_df.acc + "/" + all_features_df.pos.apply(str) + "/" + all_features_df.aa1 + "/" + all_features_df.aa2)


cur_aa["acc"] = all_features_df.acc
cur_aa["pos"] = all_features_df.pos
cur_aa["aa1"] = all_features_df.aa1
cur_aa["aa2"] = all_features_df.aa2

cur_aa = cur_aa.drop_duplicates(0)

cur_aa.set_index(0, inplace=True)

first = True
for path in tqdm(varity_paths):
    loaded_df = pd.read_csv(path)
    loaded_df[0] = loaded_df.p_vid + "/" + loaded_df.aa_pos.apply(str) + "/" + loaded_df.aa_ref + "/" + loaded_df.aa_alt
    loaded_df = loaded_df.drop_duplicates(0)
    loaded_df.set_index(0, inplace=True)
    if first:
        cur_aa = cur_aa.merge(loaded_df, how="left", left_index=True, right_index=True)
        first = False
    else:
        cur_aa.update(loaded_df)

all_features_df = all_features_df.merge(cur_aa, on=["acc", "pos", "aa1", "aa2"], how="left")

repeated_columns = ['integrated_fitCons_score_y', 'LRT_score_y',
       'phyloP30way_mammalian_y', 'phastCons30way_mammalian_y',
       'SiPhy_29way_logOdds_y']

all_features_df = all_features_df.drop(columns=repeated_columns)

all_features_df.rename(columns=dict(
    zip(list(map(lambda x: x[:-2] + "_x", repeated_columns)), list(map(lambda x: x[:-2], repeated_columns)))), inplace=True)

for field_name in constants.varity_cols:
    all_features_df.loc[:, field_name] = pd.to_numeric(all_features_df[field_name], errors="coerce")
    all_features_df[f"{field_name}_na"] = 0
    all_features_df.loc[all_features_df[field_name].isna(), f"{field_name}_na"] = 1
    all_features_df.loc[:, field_name] = all_features_df.loc[:, field_name].fillna(0)

all_features_df.to_pickle(constants.variant_features_all)