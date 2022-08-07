# pip install gtfparse


from collections import defaultdict
import csv
import glob
from gtfparse import read_gtf

gf_files = glob.glob("*.gtf")

# read positionsFile
all_loc = []
all_geno = []
dict_chr_pos_to_snp = defaultdict(list)
num_of_snps = 1
chromosoms = set()
# read chromosoms from vcf files

for file in gf_files:
    unique_set = set()
    df = read_gtf(gf_files[0])
    # filter DataFrame to gene entries on chrY
    df_genes = df[df["feature"] == "gene"]
    for row in df_genes.itertuples():
        chromosoms.add(row.seqname)
# put chromosoms in our dict
dict_of_chromosoms = {}
for index, chr in enumerate(sorted(list(chromosoms))):
    dict_of_chromosoms[chr] = "chr" + "{0:02}".format(index + 1)

for file in gf_files:
    df = read_gtf(gf_files[0])

    # filter DataFrame to gene entries on chrY
    df_genes = df[df["feature"] == "gene"]
    for row in df_genes.itertuples():
        all_loc.append((dict_of_chromosoms[row.seqname], row.start, row.end,)) if (
            dict_of_chromosoms[row.seqname],
            row.start,
            row.end,
        ) not in unique_set else None
        unique_set.add(
            (
                dict_of_chromosoms[row.seqname],
                row.start,
                row.end,
            )
        )
counter = 1
sorted_loc = sorted(all_loc, key=lambda x: (x[0], x[1]))
with open(file + "_gene_loc.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Gene", "chr", "start", "end"])
    for row in sorted_loc:
        writer.writerow(["Gene_" + str(counter), *row])
        counter += 1
