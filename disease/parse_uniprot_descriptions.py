#!/usr/bin/env python
from sys import argv
import re
script, uniprot, disease_table = argv

d_map = dict()
#Read in disease names, with keys Orphanet and OMIM IDs.
with open(disease_table) as dt:
	for line in dt:
		lines = line.strip().split("\t")
		if len(lines) > 2:
			names = lines[2].split(";")
			if (lines[0] != "Asterisk" and lines[0] != "Caret"):
				for n in names:
					d_map[n] = lines[1]

with open(uniprot) as uh:
	for line in uh:
		lines = line.strip().split("\t")
		name = lines[1].upper()
		mnemonic = lines[2]
		description = lines[3]
		if name in d_map:
			print(d_map[name] + "\t" + description)
		elif mnemonic in d_map:
			print(d_map[mnemonic] + "\t" + description)

