#  pip install scikit-allel

import csv
import pandas as pd
import allel
import glob, os

vcf_files = glob.glob("*.vcf")
# read positionsFile

for file in vcf_files:
    # read snp chrom location from file
    df = allel.vcf_to_dataframe(file, fields="*", alt_number=3)
    callset = allel.read_vcf(file, fields="*", alt_number=3)
    cromosoms = sorted(list(set(callset["variants/CHROM"])))
    dict_of_chromosoms = {}
    for index, chr in enumerate(cromosoms):
        dict_of_chromosoms[chr] = (
            "chr0" + str(index + 1) if int(index + 1) < 10 else "chr" + str(index + 1)
        )

    output = []
    for index, chrom_row in df.iterrows():
        cromosom = dict_of_chromosoms.get(chrom_row["CHROM"])
        snp = "snp0" + str(index + 1) if int(index + 1) < 10 else "snp" + str(index + 1)
        output.append(
            (
                cromosom,
                snp,
                chrom_row["POS"],
            )
        )
    with open(file + ".csv", "w") as f:
        writer = csv.writer(f)
        writer.writerow(["snp", "chromosome", "position"])
        for row in output:
            writer.writerow(row)
