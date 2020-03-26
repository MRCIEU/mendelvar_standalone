#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, genes, hgnc, out = argv

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

with open(genes) as fh:
	for line in fh:
		lines = line.strip().split("\t")
		pathway_id = lines[0].split("/")[-1]
		if pathway_id == None:
			print("NONE!")
		symbols = lines[2:]
		half = lines[1].split("; datasource")
		#Check if human
		if "9606" in half[1]:
			pass
		else:
			continue
		pathway_name = half[0].split(": ")[1]
		for s in symbols:
			my_hgnc_id = hgnc_ids.get(s.upper(), "NA")
			if my_hgnc_id == "NA":
				my_hgnc_id = old_symbols.get(s.upper(), "NA")
			if my_hgnc_id == "NA":
				print(s)
			else:
				out_h.write(my_hgnc_id + "\t" + pathway_id + "\t" + pathway_name + "\n")


out_h.close()