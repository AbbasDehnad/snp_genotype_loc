#  pip install scikit-allel

import csv
from random import sample
import pandas as pd
import allel
import glob, os

vcf_files = glob.glob("*.vcf")
# read positionsFile

for file in vcf_files:
    # read snp chrom location from file
    df = allel.vcf_to_dataframe(file, fields="*")
    # print(
    #     df.columns,
    #     df.loc[0],

    # )
    callset = allel.read_vcf(file, fields="*", alt_number=3)
    cromosoms = sorted(list(set(callset["variants/CHROM"])))
    # print(sorted(callset.keys()),allel.GenotypeArray(callset['calldata/GT']))
    dict_of_chromosoms = {}
    for index, chr in enumerate(cromosoms):
        dict_of_chromosoms[chr] = (
            "chr0" + str(index + 1) if int(index + 1) < 10 else "chr" + str(index + 1)
        )
    # TODO check snp are ordered
    snp_out = []
    unique_set = set()
    geno_out = []
    for index, chrom_row in df.iterrows():
        cromosom = dict_of_chromosoms.get(chrom_row["CHROM"])
        # GT= chrom_row["GT"]
        # print(GT)
        # add snps but dont add duplicates
        if chrom_row["is_snp"] and (cromosom, chrom_row["POS"]) not in unique_set:
            snp = (
                "snp0" + str(index + 1)
                if int(index + 1) < 10
                else "snp" + str(index + 1)
            )

            snp_out.append(
                (
                    snp,
                    cromosom,
                    chrom_row["POS"],
                )
            )
            unique_set.add((cromosom, chrom_row["POS"]))
    with open(file + "snp_loc.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["snp", "chromosome", "position"])
        # sort based on pos
        for row in sorted(list(snp_out), key=lambda x: x[2]):
            writer.writerow(row)

    with open(file + "geno_loc.csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["chr", "position", "GT"])
        geno_out = zip(
            [dict_of_chromosoms.get(chr) for chr in callset["variants/CHROM"]],
            callset["variants/POS"],
            [
                1 if gt[0][0] == 1 else 0
                for gt in allel.GenotypeArray(callset["calldata/GT"])
            ],
        )
        for row in geno_out:  # TODO check GT is correct
            writer.writerow(row)
