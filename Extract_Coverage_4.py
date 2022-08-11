#  pip install scikit-allel

from collections import defaultdict
import csv
from email.policy import default
import os
from random import sample
import allel
import glob

from pandas import DataFrame

vcf_files = glob.glob("*.vcf")
vcf_files= list(map(os.path.basename, vcf_files))
# read positionsFile
all_loc = []
all_geno = []
dict_chr_pos_to_snp = defaultdict(list)
num_of_snps = 1
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
file_num = 1
for file in vcf_files:
    print("file:", file)
    # read snp chrom location from file
    df = allel.vcf_to_dataframe(file, fields="*")
    callset = allel.read_vcf(file, fields="*", alt_number=3)
    GT = callset["calldata/GT"]
    unique_set = set()
    snp_out = []
    geno_out = []
    geno_out_1 = []
    for index, chrom_row in df.iterrows():

        cromosom = dict_of_chromosoms.get(chrom_row["CHROM"])

        snp = "snp" + str(num_of_snps)
        dict_chr_pos_to_snp[cromosom, chrom_row["POS"]].append(snp)
        num_of_snps += 1

        #  add snps but dont add duplicates
        if (cromosom, chrom_row["POS"]) not in unique_set:

            geno_out_1.append(
                (
                    (cromosom, chrom_row["POS"]),
                    1 if GT[index][0][0] == 1 else 0,
                    file_num,
                    snp,
                )
            )
            snp_out.append((snp, cromosom, chrom_row["POS"]))
        unique_set.add((cromosom, chrom_row["POS"]))

    all_loc.extend(snp_out)
    all_geno.extend(sorted(list(geno_out_1)))
    file_num += 1

# first convert it to 2d table


# sort the locations basd on chr/pos/snp
all_loc = sorted(all_loc, key=lambda x: (x[1], x[2], int(x[0][3:])))
counter_loc = 1
with open("snp_loc.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["snp", "Chromosome", "position"])
    # sort based on pos
    for row in all_loc:
        writer.writerow(["snp" + str(counter_loc), row[1], row[2]])
        counter_loc += 1

dict_snp_to_index = {snp: index for index, snp in enumerate([l[0] for l in all_loc])}
list_to_convert = [[0] * (len(vcf_files) + 1) for _ in range(num_of_snps - 1)]
for row in all_geno:
    # related snps to this chr,pos
    snp = row[3]
    snps = dict_chr_pos_to_snp[row[0]]
    list_to_convert[int(snp[3:]) - 1][0] = snp
    for snp in snps:
        col_index = row[2]
        row_index = int(snp[3:]) - 1
        list_to_convert[row_index][col_index] = row[1]


counter = 1
list_to_convert = sorted(list_to_convert, key=lambda x: dict_snp_to_index[x[0]])
with open("geno_loc.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["SNP", *vcf_files])
    for row in list_to_convert:
        writer.writerow(["snp" + str(counter), *row[1:]])
        counter += 1
