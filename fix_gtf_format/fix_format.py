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
        content_non_mt = re.sub(
            pattern=pattern_mt,
            repl="",
            string=content_new,
        )
        mt_lines = re.findall(pattern=pattern_mt, string=content_new)

    with open(f"out/converted_{file}", "w") as f:
        f.write(content_non_mt)
        for line in mt_lines:
            f.write(line)
            f.write("\n")
    with open(f"out/converted_{file}", "r") as f:
        lines= f.readlines()
    with open(f"out/converted_{file}", "w") as f:
        lines = filter(lambda x: x.strip(), lines)
        f.writelines(lines) 
  
