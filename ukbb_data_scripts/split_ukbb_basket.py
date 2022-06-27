import csv

import pandas as pd
import numpy as np

from tqdm import tqdm
from contextlib import ExitStack

import constants


basket_columns = pd.read_csv(constants.basket_path, nrows=0)
basket_phenocodes = set(list(map(lambda x: x.split("-")[0], basket_columns.columns.tolist()[1:])))

genebass_top = pd.read_csv(constants.genebass_top_hist_path)

# pick only those that we have in our list
genebass_top = genebass_top[genebass_top.Phenotype.isin(basket_phenocodes)]
genebass_phenocodes = genebass_top.Phenotype.unique().tolist()
# indices of columns that we use
use_cols = np.where(np.array([True] + list(map(lambda x: x.split("-")[0] in (genebass_phenocodes + constants.adjustment_phenocodes), basket_columns.columns.tolist()[1:]))) == True)[0]
# names of the columns
basket_cols_raw = basket_columns.columns.values[use_cols]

with open(constants.basket_path, "r") as f:
    reader = csv.reader(f)
    with ExitStack() as stack:
        files = [stack.enter_context(open(fname, "w")) for fname in 
                 [f"{constants.ukbb_extracts_output_dir}/{basket_col}.csv" for basket_col in basket_cols_raw]]
        
        # the 'total' is approximate
        for line in tqdm(reader, total=500000):
            for file, entry in zip(files, np.array(line)[use_cols]):
                file.write(entry)
                file.write('\n')
