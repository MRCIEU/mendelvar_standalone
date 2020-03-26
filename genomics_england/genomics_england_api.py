#!/usr/bin/env python
import json, requests, re
import argparse
from collections import defaultdict as dd

ap = argparse.ArgumentParser()
ap.add_argument('--genemap',required=True,type=str,help='Genemap file from OMIM')
ap.add_argument('--mimtitle',required=True,type=str,help='Processed mimTitle file with all possible disease MIM IDs')
ap.add_argument('--title',required=True,type=str,help='Original mimTitle file from OMIM')
ap.add_argument('--caret',required=True,type=str,help='Processed mimTitle file with provided mappings from disused MIM entry no to current IDs')
ap.add_argument('--genes_title',required=True,type=str,help='Processed mimTitle file with all possible gene MIM IDs')
ap.add_argument('--output',required=True,type=str,help='Output file to write parsed results')
ap.add_argument('--hgnc',required=True,type=str,help='HGNC input file')

args = ap.parse_args()

new_disease_mim = set()

#Get gene ENGS id = omim_gene id relationships
#Also gene symbol = omim_gene id relationships
gene_ids_omim = dict()
symbol_omim = dict()
seen_names = set()
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

#Read in data from OMIM title file.
#Keys: disease title, value: mim_id
mimtit_dict = {}
with open (args.title) as mimt:
    for line in mimt:
        if not line.startswith("#"):
            lines = line.strip().split("\t")
            if (lines[0].strip() != "Asterisk" and lines[0].strip() != "Caret"):
                mim = lines[1]
                title = lines[2].split(";")[0].lower().capitalize()
                mimtit_dict[title] = mim
                mimtit_dict[mim] = title

#Read in a file with old id to new id mappings (caret), so
#that we can substitute the old one with the new ID.
caret = {}
with open (args.caret) as caret_fh:
    for line in caret_fh:
        lines = line.strip().split("\t")
        caret[lines[0]] = lines[1]

disease_omim = set()
gene_omim = set()

#Read in disease MIM numbers:
with open(args.mimtitle) as tf:
    for line in tf:
        lines = line.strip().split("\t")
        disease_omim.add(lines[0])

#Read in gene MIM numbers:
with open(args.genes_title) as gtf:
    for line in gtf:
        lines = line.strip().split("\t")
        gene_omim.add(lines[0])

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

def check_moved_mim_id(omim_id):
    if omim_id in caret:
        omim_id = caret[omim_id]
    return omim_id

def check_gene_mim_exists(gene_mim):
    if gene_mim == "NA":
        return gene_mim
    elif gene_mim not in gene_omim:
        print ("Ooops, MIM gene id %s not present in OMIM!" % gene_mim) 
        return "NA"
    else:
        return gene_mim

def check_omim_id_exists(omim_id):
    if omim_id == "NA":
        return omim_id
    elif omim_id in disease_omim:
        return omim_id
    else:
        print ("Ooops, disease OMIM id from Genomics England %s not present in OMIM!" % omim_id)
        return None

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
        return my_mim

results_dict = dd(lambda: dd(dict))

success = 0
failure = 0

def parse_json(response_dict, results_dict, new_disease_mim, success, failure):
    #Mim number match
    mim_no = r".*(\d{6}).*"
    for gene in response_dict["results"]:
        omim_gene = gene["gene_data"].get("omim_gene", None)
        if omim_gene != None:
            omim_gene = omim_gene[0]
            omim_gene = check_gene_mim_exists(omim_gene)
        else:
            omim_gene = "NA"

        hgnc_symbol = gene["gene_data"].get("hgnc_symbol", "NA")
        hgnc_id =  gene["gene_data"].get("hgnc_id", "NA").replace("HGNC:", "")
        ensg_37 = gene["gene_data"]["ensembl_genes"].get("GRch37", "NA")
        if ensg_37 != "NA":
            ensg_37 = ensg_37["82"]["ensembl_id"]
        ensg_38 = gene["gene_data"]["ensembl_genes"].get("GRch38", "NA")
        if ensg_38 != "NA":
            ensg_38 = ensg_38["90"]["ensembl_id"]
        #Fill out missing gene MIM 
        if omim_gene == "NA":
            omim_gene = lookup_gene_mim(hgnc_id, ensg_38, hgnc_symbol)
            if omim_gene != "NA":
                success += 1
            else:
                failure += 1
        cf = gene.get("confidence_level", "NA")
        pheno = gene.get("phenotypes", "NA")
        for my_p in pheno:
            #Ignore phenotypes described only by MONDO or ORPHA terms (only 37 of those, and
            #mostly matching Joubert and Meckel syndrome, each of which has already >15 genes assigned)
            if (my_p.startswith("MONDO") or my_p.startswith("ORPHA")):
                continue
            my_mim_id = None
            MyResult = re.search(mim_no, my_p)
            #Match for OMIM disease ID
            if MyResult:
                my_mim_id = MyResult.group(1) 
                my_mim_id = check_moved_mim_id(my_mim_id)
                my_mim_id = check_omim_id_exists(my_mim_id)
                my_mim_name = mimtit_dict.get(my_mim_id, my_p)
            #See if we can add MIM ids to diseases with unprovided MIM disease ids.
            else:
                my_mim_name = re.sub(r"[^\w\d,\s]", "", my_p).lower().capitalize()
                if my_mim_name in mimtit_dict:
                    my_mim_id = mimtit_dict[my_mim_name]
                    new_disease_mim.add(my_mim_id)
            #Collect results into a dictionary but only if we mapped disease to OMIM id
            #And if we don't have the entry already in the dictionary with a higher confidence level:
            if my_mim_id:
                old_confidence = results_dict[my_mim_id][hgnc_id].get("status", -1)
                if int(old_confidence) > int(cf):
                    continue 
                else:
                    results_dict[my_mim_id][hgnc_id]["status"] = cf
                    results_dict[my_mim_id][hgnc_id]["omim_gene"] = omim_gene
                    results_dict[my_mim_id][hgnc_id]["hgnc_symbol"] = hgnc_symbol
                    results_dict[my_mim_id][hgnc_id]["ensg_37"] = ensg_37
                    results_dict[my_mim_id][hgnc_id]["ensg_38"] = ensg_38
                    results_dict[my_mim_id][hgnc_id]["my_mim_name"] = my_mim_name
    return results_dict, new_disease_mim, success, failure

a = requests.get("https://panelapp.genomicsengland.co.uk/api/v1/genes")
a.encoding = a.apparent_encoding

response_dict = json.loads(a.content)
results_dict, new_disease_mim, success, failure = parse_json(response_dict, results_dict, new_disease_mim, success, failure)

#Loop over all the pages of the results. Find the highest possible 
#ranking for a given gene-phenotype association

while response_dict["next"]:
    a = requests.get(response_dict["next"])
    a.encoding = a.apparent_encoding
    response_dict = json.loads(a.content)
    results_dict, new_disease_mim, success, failure = parse_json(response_dict, results_dict, new_disease_mim, success, failure)

out_fh = open(args.output, 'w')
out_fh.write("# omim_disease_id" + "\t" + "omim_disease_name" +
"\t" + "omim_gene_id" + "\t" + "hgnc_gene_name" + "\t" + "hgnc_id" +
"\t" + "ensg_gene_name_37" + "\t" + "ensg_gene_name_38" + "\n")
for my_mim_id in results_dict:
    for hgnc_id in results_dict[my_mim_id]:
        #Confidence level of at least 2 - amber
        if int(results_dict[my_mim_id][hgnc_id]["status"]) > 1:
            out_fh.write(my_mim_id + "\t" + results_dict[my_mim_id][hgnc_id]["my_mim_name"] +
            "\t" + results_dict[my_mim_id][hgnc_id]["omim_gene"] + "\t" +
            results_dict[my_mim_id][hgnc_id]["hgnc_symbol"] + "\t" +
            hgnc_id + "\t" + results_dict[my_mim_id][hgnc_id]["ensg_37"] +
            "\t" + results_dict[my_mim_id][hgnc_id]["ensg_38"] + "\n")
out_fh.close()


#Print how many diseses were successfully assigned a MIM id.
print("We found %d new disease name to OMIM disease ID mappings" % len(new_disease_mim))
print("We were able to successfully map %d of all the genes with NA gene MIM" % success)
print("We were NOT able to successfully map %d of all the genes with NA gene MIM" % failure)