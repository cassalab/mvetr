import glob
import os


vcf_files = glob.glob(constants.vcfs_dir + "/*")
# specific to our vcf creation
vcf_files = list(filter(lambda x: "ORIGINAL" not in x, vcf_files))

for vcf_filename in vcf_files:
    name = vcf_filename.rsplit("/", 1)[1].rsplit(".", 1)[0]
    output_filename = f{constants.snps_processing_dir}/{name}.mapsnps.input"
    os.system(f"./vcf2mapsnps.pl {vcf_filename} >{output_filename}")

# concatenate all of the vcfs into one big file
with open(f'{constants.snps_processing_dir}/all.mapsnps.input', 'w') as outfile:
    for fname in vcf_files:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)

# Then we run mapsnps script on the all.mapsnps.input file, e.g.: qsub -t 1-100 -cwd -b y /net/etc/mapsnps.pl -g hg38 -m -n -x1 ../mapsnps_input/all.mapsnps.input

# Finally, PolyPhen2 pipeline is executed on mysnps.mapsnps.output to extract the variant annotations
