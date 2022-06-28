import pandas as pd

import constants


# fetch top gene-phenotype pairs for traits in the data basket
def get_genebass_top():
    basket_columns = pd.read_csv(constants.basket_path, nrows=0)
    basket_phenocodes = set(list(map(lambda x: x.split("-")[0], basket_columns.columns.tolist()[1:])))

    genebass_top = pd.read_csv(constants.genebass_top_hist_path)

    # pick only those that we have in our list
    return genebass_top[genebass_top.Phenotype.isin(basket_phenocodes)]
