import pandas as pd

import constants


def parse_annotated_vep_file(input_file):
    rows_to_skip = 0
    headers = []
    with open(input_file, "r") as f:
        found = False
        line = True
        while not found and line:
            line = f.readline()
            if "#CHROM" in line and "REF" in line:
                found = True
                break
            rows_to_skip += 1
            headers.append(line)
    f.close()
    df = pd.read_csv(input_file, delimiter = "\t", skiprows = rows_to_skip)
    info_sections = headers[2].split("Ensembl VEP. Format: ")[1].replace(">\n", "").split("|")[1:]
    columns = ["Name", "Chrom", "Pos", "Ref", "Alt"] + info_sections
    rows = []

    ## 665 fields are returned by VEP with dbNSFP plugin. run using "ALL" as
    ## 
    n = 665
    for index, row in df.iterrows():
        chrom, pos, ref, alt = row["#CHROM"], int(row["POS"]), row["REF"], row["ALT"]
        if chrom not in ["X", "Y"]:
            chrom = int(chrom)
        name = "chr" + str(chrom) + ":" + str(pos) + ref + "/" + alt
        pieces =  list(map(lambda x: None if x == "" else x, row["INFO"].split("|")[1:]))
        transcripts = [pieces[i:i + n] for i in range(0, len(pieces), n)]
        canonical_fields = None
        canonical_found = 0
        for c in transcripts:
            d = {}
            for field, value in zip(info_sections, c):
                d[field] = value

            if d["CANONICAL"] and "YES" in d["CANONICAL"]:
                canonical_fields = c
                canonical_found = 1
                break

        if canonical_fields is None:
            canonical_fields = transcripts[0]
        assert len(pieces) % n == 0
        rows.append([name, chrom, pos, ref, alt] + canonical_fields + [canonical_found])

    variant_df = pd.DataFrame(rows, columns = columns + ["canonical_found"])
    return variant_df

vcf_output = parse_annotated_vep_file(constants.vcf_output_path)

vcf_output = vcf_output[vcf_output.BIOTYPE == "protein_coding"]

vcf_output = vcf_output[vcf_output.Consequence.apply(lambda x: "missense_variant" in x)]

vcf_output = vcf_output[["Name"] + constants.dbnsfp_continuous_cols + constants.dbnsfp_categorical_cols]

vcf_output.rename(columns={"Name": "var_id"}, inplace=True)

for field_name in constants.dbnsfp_continuous_cols:
    vcf_output.loc[:, field_name] = pd.to_numeric(vcf_output[field_name], errors="coerce")
    vcf_output[f"{field_name}_na"] = 0
    vcf_output.loc[vcf_output[field_name].isna(), f"{field_name}_na"] = 1
    vcf_output.loc[:, field_name] = vcf_output.loc[:, field_name].fillna(0)

for col_name in constants.dbnsfp_categorical_cols:
    vcf_output[col_name] = vcf_output[col_name].astype(str).apply(str.strip).astype("category").cat.codes

# merge with the other features

all_features_df = pd.read_pickle(constants.variant_features_all)

vcf_output = vcf_output.merge(all_features_df, on="var_id", how="left")

vcf_output.to_pickle(constants.variant_features_all)
