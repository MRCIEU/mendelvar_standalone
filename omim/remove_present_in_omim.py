#!/usr/bin/env python
import argparse, json
from collections import defaultdict as dd
#The script compares the gene-disease relationships provided by the file, 
#and removes those which duplicate the ones present in the OMIM genemap table.

ap = argparse.ArgumentParser()
ap.add_argument('--omim',required=True,type=str,help='Output from parse_genemap_relationships.py script (OMIM)')
ap.add_argument('--other',required=True,type=str,help='Another database file for comparison')
ap.add_argument('--omim_gene',required=True,type=str,help='Column name for gene id in OMIM file')
ap.add_argument('--omim_disease',required=True,type=str,help='Column name for disease id in OMIM file')
ap.add_argument('--other_gene',required=True,type=str,help='Column name for gene id in the other file')
ap.add_argument('--other_disease',required=True,type=str,help='Column name for disease id in the other file')
ap.add_argument('--pheno_series',required=True,type=str,help='A json dictionary containing mappings between disease mim IDs ant their phenotypic series number')
args = ap.parse_args()


omim_dict = dd(set)
with open (args.omim) as omim_fh:
	header = omim_fh.readline().strip().split("\t")
	for line in omim_fh:
		lines = line.strip().split("\t")
		#Dictionary indexed by gene id, containing a set of disease identifiers as values
		omim_dict[lines[header.index(args.omim_gene)]].add(lines[header.index(args.omim_disease)])

with open (args.pheno_series, 'r') as ps_h:
	ps = json.load(ps_h)


with open (args.other) as other_fh:
	header = other_fh.readline().strip().split("\t")
	print ("\t".join(header))
	for line in other_fh:
		lines = line.strip().split("\t")
		#Check if our gene is present in the OMIM data
		if lines[header.index(args.other_gene)] in omim_dict:
			#Check if our disease ID is present in OMIM for our gene ID
			if lines[header.index(args.other_disease)] in omim_dict[lines[header.index(args.other_gene)]]:
				pass
			#If not, test if perhaps the disease is part of phenotypic series, and the gene-disease link is still present
			#in OMIM under a disease ID in the same phenotypic series.
			elif lines[header.index(args.other_disease)] in ps:
				to_ignore = 0
				for alt_ids in ps[lines[header.index(args.other_disease)]]:
					if str(alt_ids) in omim_dict[lines[header.index(args.other_gene)]]:
						to_ignore = 1
				#If gene-disease association not present in OMIM, print out
				if to_ignore == 0:
					print(line.rstrip())
			#If gene-disease association not present in OMIM, print out
			else:
				print(line.rstrip())
		#If gene absent in OMIM data, print out
		else:
			print(line.rstrip())
