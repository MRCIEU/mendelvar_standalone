#!/usr/bin/env python
from sys import argv

#Parse the GO annotation file to bring it in line with that of DO and HPO GAFs.
script, gaf, gpi, out = argv


hgnc_dict = dict()
with open(gpi) as gpi_h:
	for line in gpi_h:
		if line.startswith("!"):
			pass
		else:
			lines = line.strip().split("\t")
			if lines[8].startswith("HGNC:"):
				hgnc_id = lines[8].replace("HGNC:", "")
				hgnc_dict[lines[1]] = hgnc_id

out_h = open(out, 'w')
		
with open(gaf) as gaf_h:
	for line in gaf_h:
		if line.startswith("!"):
			out_h.write(line)
		else:
			lines = line.strip().split("\t")
			if (lines[3].startswith("NOT") or lines[6] == "IEA"):
				pass
			else:
				if lines[1] in hgnc_dict:
					lines[0] = "HGNC"
					lines[1] = hgnc_dict[lines[1]]
					lines = lines[0:15]
					out_h.write("\t".join(lines) + "\n")

out_h.close()