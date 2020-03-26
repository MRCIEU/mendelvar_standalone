#!/usr/bin/env python
import glob, os
from sys import argv
from collections import defaultdict as dd

script, genes, hgnc, out = argv

to_glob = genes + "/*bed"

hgnc_ids = dict()
old_symbols = dict()
with open(hgnc, 'rb') as hgnc_h:
    header = hgnc_h.readline().decode('utf8').strip().split("\t")
    for line in hgnc_h:
        lines = line.decode('utf8').strip().split("\t")
        hg_id = lines[0].replace("HGNC:", "")
        hg_symbol = lines[1]
        hgnc_ids[hg_symbol.upper()] = hg_id
        alias = lines[8].replace('"', '').strip()
        prev_symbol = lines[10].replace('"', '').strip()
        if len(alias) > 0:
        	aliases = alias.split("|")
        	for a in aliases:
        		old_symbols[a.upper()] = hg_id
        if len(prev_symbol) > 0:
        	prev_symbols = prev_symbol.split("|")
        	for p in prev_symbols:
        		old_symbols[p.upper()] = hg_id

out_h = open(out, "w")

counter = 0
for filename in glob.iglob(to_glob):
	with open(filename) as fh:
		counter += 1
		for line in fh:
			lines = line.strip().split("\t")
			my_file = os.path.basename(filename)
			category = os.path.splitext(my_file)[0]
			gene_set_id = "Freund_" + str(counter)
			my_hgnc_id = hgnc_ids.get(lines[3].upper(), "NA")
			if my_hgnc_id == "NA":
				my_hgnc_id = old_symbols.get(lines[3].upper(), "NA")
			if my_hgnc_id == "NA":
				continue
			else:
				out_h.write(my_hgnc_id + "\t" + gene_set_id + "\t" + category + "\n")


out_h.close()