# Missense Variant Trait Effect predictor (MVTE)
Inspired by MVP, we develop a deep learning network to predict missense variant effects on traits in humans. We use exome sequencing and trait data from UKBB.

## Steps
### 0. Update the paths
Update the paths in the `ukbb_data_scripts/constants.py` file.

### 1. Phenotype dataset preprocessing

#### a) Identify gene-phenotype associations
We limit our study to the top gene-phenotype associations identified in [Systematic single-variant and gene-based association testing of thousands of phenotypes in 426,370 UK Biobank exomes](https://www.medrxiv.org/content/10.1101/2021.06.19.21259117v4.full-text). We take the results of the SKAT-O gene-based burden tests with p-value threshold of 10^{-5}. They can be downloaded using the web-GUI at [https://app.genebass.org/gene/ENSG00000175445/phenotype/continuous-30760-both_sexes--irnt?burdenSet=missense%7CLC&phewasOpts=1&resultIndex=top-associations&resultLayout=full](https://app.genebass.org/gene/ENSG00000175445/phenotype/continuous-30760-both_sexes--irnt?burdenSet=missense%7CLC&phewasOpts=1&resultIndex=top-associations&resultLayout=full). Update the `genebass_top_hist_path` in `ukbb_data_scripts/constants.py`.

#### b) Extract the relevant phenotype annotations
Execute `ubb_data_scripts/split_ukbb_basket.py`

The script extracts the phenothype annotations for the phenotypes identified in the previous step that appear in the data basket and two more for trait adjustment: [Age when attended assessment centre](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=21003), and [Sex](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=31).

This step reduces the file sizes for ease of computing.

Similarly extract [10 genetic principal components](https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22009) to a csv with "eid" column that corresponds to the patient ids. Update the  `genetic_pcs_path` in `ukbb_data_scripts/constants.py`.

#### c) Generate patient_id-phenotype dataframes
To create a "statin" phenotype that indicates whether a subject takes any of the listed statins, download the medication coding https://biobank.ndph.ox.ac.uk/showcase/coding.cgi?id=4 and update the `medication_coding_path` in `ukbb_data_scripts/constants.py`.

Execute `ukbb_data_scripts/create_eid_trait_pickles.py`

The script parses the csvs generated in step 1. b) to create continuous and binary trait dataframes.

#### d) Adjust for covariates
Execute `ukbb_data_scripts/adjust_for_covariates.py`

Adjusts for sex, age, 10PCs, age^2, sex*age, sex*age^2 using linear regression for the continuous traits. Lipid metabolites ("Cholesterol", "LDL direct", "Apolipoprotein B", "HDL cholesterol", "Apolipoprotein A") are also adjusted for the intake of statins. Binary traits are left untreated. In addition, it creates gene-phenotype association dataframe.

### 2. Variant dataset preprocessing

#### a) Filter missense variants, extract pp2 annotations
We write the UKBB variants to vcf format and store the path to the directory in `ukbb_data_scripts/constants.py` in `vcfs_dir` variable. Note that these scripts use [PolyPhen2 tools](http://genetics.bwh.harvard.edu/pph2/dokuwiki/downloads) to extract variant annotations listed [here](http://genetics.bwh.harvard.edu/pph2/dokuwiki/appendix_a). Our pipeline is described in `ukbb_data_scripts/polyphen2_pipeline.py`. The PP2 annotations are stored in `pp2_output_path`,  found in `ukbb_data_scripts/constants.py`.

#### b) Add EVE and PP2 scores
Download [EVE scores](https://evemodel.org/download/bulk) and update the `eve_scores_dir` in `ukbb_data_scripts/constants.py` to point to the `variant_files` directory. `vcf_input_path` and `vcf_output_path` should be set too, that's where vcf input/output to/from [vep](http://uswest.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#dbnsfp) will be. Execute `ukbb_data_scripts/extract_features.py`

<!-- Same for [VARITY](http://varity.varianteffect.org/) and the `varity_scores_path` variable. -->

#### c) Extract dbNSFP annotations
Execute vep on the input in VCF format stored in `vcf_input_path` and store the output in `vcf_output_path`, e.g. by running `vep --cache --offline --dir /net/data/vep --assembly GRCh38 --plugin dbNSFP,/net/home/tianyu/dbNSFP4.2a_grch38.gz,ALL --vcf --no_check_variants_order --canonical --no_stats -I <vcf_input_path> -o <vcf_output_path>`.
To parse and add the dbNSFP annotations, execute `ukbb_data_scripts/add_dbnsfp_features.py`. In this step we also limit the variants to those idenitified as missense by dbNSFP.

#### d) Add VARITY annotations
Execute `ukbb_data_scripts/add_varity_features.py` to download [VARITY features](http://varity.varianteffect.org/) and add them to the dataframe with variant features. Note that this requires selenium and chrome web driver.

#### e) Get variant statistics
In this step we generate a dataset of variant descriptors. We save 

