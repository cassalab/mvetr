# Missense Variant Trait Effect predictor (MVTE)
Inspired by MVP, we develop a deep learning network to predict missense variant effects on traits in humans. We use exome sequencing and trait data from UKBB.

## Steps
### 0. Update the paths
Change the paths in the ukbb_data_scripts/constants.py file.

### 1. Identify gene-phenotype associations
We limit our study to the top gene-phenotype associations identified in [Systematic single-variant and gene-based association testing of thousands of phenotypes in 426,370 UK Biobank exomes](https://www.medrxiv.org/content/10.1101/2021.06.19.21259117v4.full-text). We take the results of the SKAT-O gene-based burden tests with p-value threshold of 10^{-5}. They can be downloaded using the web-GUI at [https://app.genebass.org/gene/ENSG00000175445/phenotype/continuous-30760-both_sexes--irnt?burdenSet=missense%7CLC&phewasOpts=1&resultIndex=top-associations&resultLayout=full](https://app.genebass.org/gene/ENSG00000175445/phenotype/continuous-30760-both_sexes--irnt?burdenSet=missense%7CLC&phewasOpts=1&resultIndex=top-associations&resultLayout=full). Update the genebass_top_hist_path in ukbb_data_scripts/constants.py.

### 2. Extract the relevant phenotype annotations
Execute `ubb_data_scripts/split_ukbb_basket.py`

The script extracts the phenothype annotations for the phenotypes identified in the previous step that appear in the data basket and two more for trait adjustment: [Age when attended assessment centre](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003), and [Sex](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31).

This step reduces the file sizes for ease of computing.

### 3. Generate patient_id-phenotype numpy arrays

