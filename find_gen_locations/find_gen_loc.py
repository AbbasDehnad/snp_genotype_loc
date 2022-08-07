#  pip install scikit-allel
import re
import csv
import glob


fna_files = glob.glob("*.fna")
print(fna_files)
from Bio import SeqIO


for file in fna_files:
    dict_chr = {}
    output = []
    counter = 1
    num_of_cromosoms = 1
    existing_snps = set()
    with open(file, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            name = record.name.split("_")[0].split("|")[1]
            chr = dict_chr.get(name, None)
            if chr is None:
                dict_chr[name] = "chr" + str(num_of_cromosoms)
                num_of_cromosoms += 1
            desc = record.description.split(" ")[
                5
            ]  # see if "location is in the 5th position"
            matches = re.findall("[0-9.]+", desc)
            if (
                not matches or ".." not in matches[0]
            ):  # should be in this format (6377..7226)
                descs = record.description.split(" ")
                desc = [s for s in descs if "location" in s][
                    0
                ]  # get the location description from uncertain places
                matches = re.findall("[0-9.]+", desc)
            if len(matches) == 1:
                try:
                    start = int(matches[0].split("..")[0])
                    end = matches[0].split("..")[1]
                except:
                    print("failed", matches, record.description)
                    continue
            else:
                start = int(matches[0].split("..")[0])
                temp_end = matches[-1].split("..")
                end = int(temp_end[1]) if len(temp_end) == 2 else int(temp_end[0])
            if (dict_chr[name], start, end) not in existing_snps:
                existing_snps.add((dict_chr[name], start, end))
                output.append([dict_chr[name], start, end])

    # sort output by chr and start
    output = sorted(output, key=lambda x: (x[0], int(x[1])))
    # write list output to csv file
    with open(
        file + ".csv",
        "w",
        newline="",
    ) as f:
        writer = csv.writer(f)
        writer.writerow(["Gene", "chr", "start", "end"])
        counter = 1
        for row in output:
            snp = "Gene_" + str(counter)
            writer.writerow([snp, *row])
            counter += 1
