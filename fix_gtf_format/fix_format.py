import re

# pip install gtfparse


from collections import defaultdict
import csv
import glob
from gtfparse import read_gtf

gf_files = glob.glob("*.gtf")
patern1 = "DBxref.*;\sID"
patern2 = 'Ontology_term.*]";'
pattern_mt = "mtDNA.*;"
pattern_not_mt = "^[^mt].*"
note_pattern = 'note .*original'
product_pattern = '" " .*;'
chr1="CP031385.1"
chr2="CP031386.1"
chr3="CP031387.1"
chr4="CP031388.1"
chr5="CP031389.1"
chr6="CP031390.1"
chr7="CP031391.1"
mtDNA="CP031392.1"
chrom_list = [chr1, chr2, chr3, chr4, chr5, chr6, chr7, mtDNA]

for file in gf_files:
    with open(file, "r") as f:
        content = f.read()
        content_new = re.sub(
            pattern=patern1,
            repl="ID",
            string=content,
        )
        content_new = re.sub(
            pattern=patern2,
            repl="",
            string=content_new,
        )
        content_new = re.sub(
            pattern=note_pattern,
            repl="original",
            string=content_new,
        )
        content_new = re.sub(
            pattern=product_pattern,
            repl='";',
            string=content_new,
        )
        content_new = re.sub(
            pattern=note_pattern,
            repl="original",
            string=content_new,
        )
        content_non_mt = re.sub(
            pattern=pattern_mt,
            repl="",
            string=content_new,
        )
        mt_lines = re.findall(pattern=pattern_mt, string=content_new)

    with open(f"converted_{file}", "w") as f:
        f.write(content_non_mt)
        for line in mt_lines:
            f.write(line)
            f.write("\n")
    with open(f"converted_{file}", "r") as f:
        content_new = f.read()
        for chr,repl in zip(["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","mtDNA" ],chrom_list):
            content_new = re.sub(
                pattern=chr,
                repl=repl,
                string=content_new,
            )
    with open(f"converted_{file}", "w") as f:
        f.write(content_new)

    with open(f"converted_{file}", "r") as f:
        lines= f.readlines()
    with open(f"converted_{file}", "w") as f:
        lines = filter(lambda x: x.strip(), lines)
        f.writelines(lines) 
  
