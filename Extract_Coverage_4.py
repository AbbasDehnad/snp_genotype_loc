#  pip install scikit-allel

import csv
from random import sample
import allel
import glob

vcf_files = glob.glob("*.vcf")
# read positionsFile
all_loc = []
all_geno = []
chromosoms = set()
# read chromosoms from vcf files
for file in vcf_files:
    callset = allel.read_vcf(file, fields="*", alt_number=3)
    cromosom = set(callset["variants/CHROM"])
    chromosoms = chromosoms.union(cromosom)
# put chromosoms in our dict
dict_of_chromosoms = {}
for index, chr in enumerate(sorted(list(chromosoms))):
    dict_of_chromosoms[chr] = "chr" + "{0:02}".format(index + 1)

print(dict_of_chromosoms)

for file in vcf_files:
    # read snp chrom location from file
    df = allel.vcf_to_dataframe(file, fields="*")
    callset = allel.read_vcf(file, fields="*", alt_number=3)

    unique_set: set[str, str] = set()
    snp_out = []
    geno_out = []
    for index, chrom_row in df.iterrows():
        cromosom = dict_of_chromosoms.get(chrom_row["CHROM"])

        # add snps but dont add duplicates
        if chrom_row["is_snp"] and (cromosom, chrom_row["POS"]) not in unique_set:
            snp = "snp" + "{0:03}".format(index + 1)
            snp_out.append((snp, cromosom, chrom_row["POS"], file))
        unique_set.add((cromosom, chrom_row["POS"]))
        geno_out = zip(
            [dict_of_chromosoms.get(chr) for chr in callset["variants/CHROM"]],
            callset["variants/POS"],
            [
                1 if gt[0][0] == 1 else 0
                for gt in allel.GenotypeArray(callset["calldata/GT"])
            ],
            [file for _ in callset["calldata/GT"]],
        )
    all_loc.extend(sorted(list(snp_out), key=lambda x: (x[1], x[2])))
    all_geno.extend(sorted(list(geno_out), key=lambda x: (x[0], x[1])))


with open("snp_loc.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["snp", "Chromosome", "position", "file_name"])
    # sort based on pos
    for row in all_loc:
        writer.writerow(row)

with open("geno_loc.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Chromosome", "position", "GT", "file_name"])

    for row in all_geno:  # TODO check GT is correct
        writer.writerow(row)
