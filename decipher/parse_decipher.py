#!/usr/bin/env python
import csv, re, argparse

ap = argparse.ArgumentParser()
ap.add_argument('--delimit',required=True,type=str,help='Delimiter used in the file')
ap.add_argument('--input',required=True,type=str,help='Input file')

args = ap.parse_args()

hpo_re = r"HP\:\d+"

if args.delimit == "csv":
	delimit = ","
else:
	delimit = "\t"

with open(args.input) as csvfile:
	csv_reader = csv.reader(csvfile, delimiter=delimit, quotechar='"')
	header = next(csv_reader) 
	sub_header = ["hgnc_gene_name", "omim_gene_id", "disease_name", "omim_disease_id", "hpo", "organ_specificity_list", "hgnc_id"]

	print('\t'.join(sub_header))
	for row in csv_reader:
		if (row[4] == "possible"):
			pass
		else:
			all_matches_hpo = list()
			all_matches_hpo.extend(re.findall(hpo_re, row[7]))
			#Unique HPO terms
			all_matches_hpo = [a.replace("HP:", "") for a in all_matches_hpo]
			all_matches_hpo = set(all_matches_hpo)	
			row[7] = ";".join(all_matches_hpo)
			sub_row = [row[0], row[1], row[2], row[3], row[7], row[8], row[12]]
			print('\t'.join(sub_row))

