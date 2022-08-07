#  pip install scikit-allel
# pip install bcbio-gff
from collections import defaultdict
import csv
import glob
from BCBio import GFF


gf_files = glob.glob("*.gff")
# read positionsFile
all_loc = []
all_geno = []
dict_chr_pos_to_snp = defaultdict(list)
num_of_snps = 1
chromosoms = set()
# read chromosoms from vcf files
for file in gf_files:
    in_handle = open(file)
    for req in GFF.parse(in_handle):
        chromosoms.add(req.id)
    in_handle.close()
# put chromosoms in our dict
dict_of_chromosoms = {}
for index, chr in enumerate(sorted(list(chromosoms))):
    dict_of_chromosoms[chr] = "chr" + "{0:02}".format(index + 1)

print(dict_of_chromosoms)
file_num = 1
for file in gf_files:
    print("file:", file)
    # read snp chrom location from file
    in_handle = open(file)

    unique_set: set[str, str] = set()
    snp_out = []
    geno_out = []
    geno_out_1 = []
    for req in GFF.parse(in_handle):
        if file.endswith("1.gTf"):
            print(req.id, req.start, req.end)
        for feature in req.features:
            if feature.type == "gene":
                geno_out.append(
                    (
                        dict_of_chromosoms[req.id],
                        feature.location.start,
                        feature.location.end,
                    )
                ) if (
                    dict_of_chromosoms[req.id],
                    feature.location.start,
                    feature.location.end,
                ) not in unique_set else None
                unique_set.add(
                    (
                        dict_of_chromosoms[req.id],
                        feature.location.start,
                        feature.location.end,
                    )
                )

    counter = 1
    list_to_convert = sorted(geno_out, key=lambda x: (x[0], int(x[1])))
    with open(file + "_loc.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", "chr", "start", "end"])
        for row in list_to_convert:
            writer.writerow(["Gene_" + str(counter), *row])
            counter += 1
    in_handle.close()
