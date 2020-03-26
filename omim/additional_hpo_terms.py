#!/usr/bin/env python
import argparse, json
from collections import defaultdict as dd

#Compare two files (from OMIM and some other data source) and HPO, DO etc. terms matched to a given disease MIM in each file.
#Write lines from the other data source file just with HPO, DO etc. terms not present in the OMIM file, any other fields present.
#Print the total number of additional HPO terms.

ap = argparse.ArgumentParser()
ap.add_argument('--omim',required=True,type=str,help='Output from parse_genemap_relationships.py script (OMIM)')
ap.add_argument('--other',required=True,type=str,help='Another database file for comparison')
ap.add_argument('--omim_hpo',required=True,type=str,help='Column name for HPO ids in OMIM file (Ids separated by semi-colons)')
ap.add_argument('--omim_disease',required=True,type=str,help='Column name for disease id in OMIM file')
ap.add_argument('--other_hpo',required=True,type=str,help='Column name for HPO ids in the other file (Ids separated by semi-colons)')
ap.add_argument('--other_disease',required=True,type=str,help='Column name for disease id in the other file')
ap.add_argument('--pheno_series',required=True,type=str,help='A json dictionary containing mappings between disease mim IDs ant their phenotypic series number')
ap.add_argument('--output',required=True,type=str,help='A name for the output file')
args = ap.parse_args()

omim_dict = dd(set)
with open (args.omim) as omim_fh:
	header = omim_fh.readline().strip().split("\t")
	for line in omim_fh:
		lines = line.strip().split("\t")
		#Dictionary indexed by disease omim id, containing a set of hpo terms as values
		all_hpo = lines[header.index(args.omim_hpo)].split(";")
		omim_dict[lines[header.index(args.omim_disease)]].update(set(all_hpo))

with open (args.pheno_series, 'r') as ps_h:
	ps = json.load(ps_h)

output_fh = open(args.output, 'w')

#Count the number of added HPO terms
hpo_counter = 0
with open (args.other) as other_fh:
	header = other_fh.readline().strip().split("\t")
	output_fh.write(("\t".join(header)) + "\n")
	for line in other_fh:
		lines = line.strip().split("\t")
		#Check if our disease is present in the OMIM data
		if lines[header.index(args.other_disease)] in omim_dict:
			#Check if our HPO terms are present in the OMIM data.
			omim_hpo = omim_dict[lines[header.index(args.other_disease)]]
			other_hpo = set(lines[header.index(args.other_hpo)].split(";"))
			unique_other = other_hpo.difference(omim_hpo) 
			#Remove any "NA" terms
			unique_other.discard("NA")
			if len(unique_other) > 0:
				hpo_counter = hpo_counter + len(unique_other)
				lines[header.index(args.other_hpo)] = ";".join(unique_other)
				output_fh.write(("\t".join(lines)) + "\n")
		else:
			output_fh.write(("\t".join(lines)) + "\n")

print("We found %d of new terms in the other file" % hpo_counter)

output_fh.close()