#!/usr/bin/env python
import argparse, csv, re
#Due to multiple issues with Orphanet data, going to evaluate the correctness of Orphanet data:
#check if MIM gene IDS given are correct for the name (check_gene_mim_matches, check passed)
#check if some missing MIM gene IDS can be supplanted by looking up HGCN names in HGCN and cross-referencing
#with ENSG ids which can then be used in OMIM.
#check if all the MIM gene IDS (check_gene_mim_exist, check passed) and disease IDS exist in OMIM (check_omim_id_exists, check failed).

ap = argparse.ArgumentParser()
ap.add_argument('--genemap',required=True,type=str,help='Genemap file from OMIM')
ap.add_argument('--title',required=True,type=str,help='Processed mimTitle file with all possible disease MIM IDs')
ap.add_argument('--genes_title',required=True,type=str,help='Processed mimTitle file with all possible gene MIM IDs')
ap.add_argument('--caret',required=True,type=str,help='Processed mimTitle file with provided mappings from disused MIM entry no to current IDs')
ap.add_argument('--orphanet',required=True,type=str,help='Orphanet input file')
ap.add_argument('--hgnc',required=True,type=str,help='HGNC input file')
ap.add_argument('--output',required=True,type=str,help='Name for output file with corrected Orphanet table, which includes added gene MIM ids')

args = ap.parse_args()

gene_ids_omim = dict()
symbol_omim = dict()
seen_names = set()
disease_omim = set()
gene_omim = set()
caret = dict()
#Get gene ENGS id = omim_gene id relationships
#Also gene symbol = omim_gene id relationships
with open(args.genemap) as gf:
	for line in gf:
		if line.startswith("#"):
			next
		else:
			lines = line.strip().split("\t")
			if (len(lines) > 10 and lines[10] != "" and lines[5] != ""):
				if lines[10] in seen_names:
					print("Oops, the name %s is not unique" % lines[10])
				else:
					seen_names.add(lines[10])
					gene_ids_omim[lines[10]] = lines[5]
			if (len(lines) > 8 and lines[8] != "" and lines[5] != ""):
				symbol_omim[lines[8]] = lines[5]

#Read in disease MIM numbers:
with open(args.title) as tf:
	for line in tf:
		lines = line.strip().split("\t")
		disease_omim.add(lines[0])

#Read in gene MIM numbers:
with open(args.genes_title) as gtf:
	for line in gtf:
		lines = line.strip().split("\t")
		gene_omim.add(lines[0])

#Read in moved MIM numbers with caret:
with open(args.caret) as cf:
	for line in cf:
		lines = line.strip().split("\t")
		caret[lines[0]] = lines[1]

gene_ids_hgnc = dict()
gene_ids_symbol = dict()
symbol_to_hgnc = dict()
symbol_to_ensg = dict()
#Get hgnc id = ensg id relationship
#Also hgnc id = gene symbol relationships
#Also gene symbol = hgnc id, ensg_id relationship
with open(args.hgnc) as hf:
	header = hf.readline()
	for line in hf:
		lines = line.strip().split("\t")
		hg_id = lines[0].replace("HGNC:", "")
		if (len(lines) > 19 and lines[19] != "" and lines[0] != ""):	
			gene_ids_hgnc[hg_id] = lines[19]
			symbol_to_ensg[lines[1]] = lines[19]
		gene_ids_symbol[hg_id] = lines[1]
		symbol_to_hgnc[lines[1]] = hg_id
gene_ids_hgnc["10604"] = "ENSG00000284194"
gene_ids_hgnc["25244"] = "ENSG00000284862"

#Read in Orphanet data.
def check_gene_mim_exist(gene_mim):
	if gene_mim == "NA":
		return gene_mim
	elif gene_mim not in gene_omim:
		print ("Ooops, MIM gene id %s not present in OMIM!" % gene_mim) #Check passed
		return "NA"
	else:
		return gene_mim
def lookup_gene_mim(hgcn_id, my_eng, symbol):
		#If no ENSG present in input
		if my_eng == "NA":
			my_eng = gene_ids_hgnc.get(hgcn_id, "NA")
			my_mim = gene_ids_omim.get(my_eng, "NA")
		else:
			my_mim = gene_ids_omim.get(my_eng, "NA")
		#If still empty after matching by ENSG, match by symbol
		if my_mim == "NA":
			my_hgnc = symbol_to_hgnc.get(symbol, "NA")
			if my_hgnc != "NA":
				my_eng = gene_ids_hgnc.get(my_hgnc, "NA")
				my_mim = gene_ids_omim.get(my_eng, "NA")
		return my_mim, my_eng


def check_omim_id_exists(omim_id):
	if omim_id == "NA":
		return omim_id
	elif omim_id in disease_omim:
		return omim_id
	else:
		print ("Ooops, disease OMIM id from Orphanet %s not present in OMIM!" % omim_id)
		return "NA"

def check_moved_mim_id(omim_id):
	if omim_id in caret:
		omim_id = caret[omim_id]
	return omim_id


success = 0
failure = 0
out = open (args.output, 'w')
with open(args.orphanet) as csvfile:
	csv_reader = csv.reader(csvfile, delimiter='\t', quotechar='"')
	header = next(csv_reader) 
	header[1] = "omim_disease_id"
	header[4] = "hgnc_gene_name"
	header[5] = "ensg_gene_name"
	header[6] = "hgnc_id"
	header[7] = "omim_gene_id"
	header[9] = "hpo"

	out.write("\t".join(header) + "\n")
	for row in csv_reader:
		#Check if disease or gene MIM can be mapped to one which has since been replaced.
		row[1] = check_moved_mim_id(row[1])
		row[7] = check_moved_mim_id(row[7])
		row[7] = check_gene_mim_exist(row[7]) #Check if gene MIM exists in OMIM
		row[1] = check_omim_id_exists(row[1]) #Check if disease MIM exists in OMIM
		#Fill out missing gene MIM and ENSG data
		if row[7] == "NA":
			row[7], row[5] = lookup_gene_mim(row[6], row[5], row[4])
			if row[7] != "NA":
				success += 1
			else:
				failure += 1
		out.write("\t".join(row) + "\n")

out.close()
print("We were able to successfully map %d of all the genes with NA gene MIM" % success)
print("We were NOT able to successfully map %d of all the genes with NA gene MIM" % failure)