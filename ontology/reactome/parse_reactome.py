#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, genes, hgnc, uniprot, out = argv

hgnc_ids = dict()
with open(hgnc, 'rb') as hgnc_h:
    header = hgnc_h.readline().decode('utf8').strip().split("\t")
    for line in hgnc_h:
        lines = line.decode('utf8').strip().split("\t")
        hg_id = lines[0].replace("HGNC:", "")
        hg_symbol = lines[1]
        if (len(lines) > 19 and lines[19] != "" and lines[0] != ""):    
            ensg_id = lines[19]
            hgnc_ids[ensg_id] = hg_id
        if (len(lines) > 25 and lines[25] != "" and lines[0] != ""): 
        	uni_id = lines[25]
        	hgnc_ids[uni_id] = hg_id
        #Try mapping by omim ID or Uniprot ID
        if (len(lines) > 31 and lines[31] != "" and lines[0] != ""):  
        	mim_id = lines[31]
        	hgnc_ids[mim_id] = hg_id

#In same cases, reactome gives ENST rather than ENSG. Read in ENST to ENSG mappings from UniProt.
#Reactome file contains alternative gene alleles. See if we have one representative of the gene.
mim_dict = dict()
uniprot_ids_dict = dict()
uniprot_dict_tr = dd(set)
with open(uniprot, 'r') as up_h:
  for line in up_h:
  	lines = line.strip().split("\t")
  	if len(lines) > 19:
  		ensg = lines[18]
  		enst = lines[19]
  		mim = lines[13]
  		ensgs = [e.strip() for e in ensg.split(";")]
  		ensts = [e.strip() for e in enst.split(";")]
  		for e in ensgs:
  			uniprot_ids_dict[e] = lines[0]
  			if len(mim) > 5:
  				mim_dict[e] = mim
  			for t in ensts:
  				uniprot_dict_tr[t].add(e)
  				uniprot_ids_dict[t] = lines[0]
  				if len(mim) > 5:
  					mim_dict[t] = mim


def get_hgnc(ens_id):
	if ens_id in hgnc_ids:
		return hgnc_ids[ens_id]
	if ens_id in uniprot_dict_tr:
		gene_ensg = uniprot_dict_tr[ens_id]
		for g in gene_ensg:
			if g in hgnc_ids:
				return hgnc_ids[g]
			if g in mim_dict:
				my_mim = mim_dict[g]
				if my_mim in hgnc_ids:
					return hgnc_ids[my_mim]
			if g in uniprot_ids_dict:
				my_uni = uniprot_ids_dict[g]
				if my_uni in hgnc_ids:
					return hgnc_ids[my_uni]
	if ens_id in mim_dict:
		my_mim = mim_dict[ens_id]
		if my_mim in hgnc_ids:
			return hgnc_ids[my_mim]	
	if ens_id in uniprot_ids_dict:
		my_uni = uniprot_ids_dict[ens_id]
		if my_uni in hgnc_ids:
			return hgnc_ids[my_uni]
	return "NA"


out_h = open(out, "w")
present = open("ens_present.txt", "w")

all_present = set()
with open(genes) as fh:
	for line in fh:
		lines = line.strip().split("\t")
		#Check if human
		if lines[5] == "Homo sapiens" and lines[4] != "IEA":
			pass
		else:
			continue
		ens_id = lines[0]
		pathway_id = lines[1]
		pathway_name = lines[3]
		hg_id = get_hgnc(ens_id)
		if hg_id == "NA":
			print(lines[0])
		else:
			out_h.write(hg_id + "\t" + pathway_id + "\t" + pathway_name + "\n")
			all_present.add(ens_id)

for a in all_present:
	present.write(a + "\n")

out_h.close()
present.close()

