#!/usr/bin/env python
from sys import argv
import re
from collections import defaultdict as dd

script, omim, orphanet, decipher, ge, hgnc, out = argv

#Store hgnc id matches to ENSG
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
            hgnc_ids[hg_id] = ensg_id
            hgnc_ids[hg_symbol] = hg_id
        else:
            hgnc_ids[hg_symbol] = hg_id

omim_data =  dd(lambda: dd(dict))
#keys: gene_hgnc_id, disease_omim_id (or disease name if absent)
omim_data_by_disease = dd(lambda: dd(dict))
omim_names = dd()
#Key: orphanet disease id, value: omim disease id:
orphanet_mappings = dd(lambda: dd(int))
with open(omim, 'rb') as omim_h:
    header = omim_h.readline().decode('utf8').strip().split("\t")
    for line in omim_h:
        lines = line.decode('utf8').strip().split("\t")
        #Get hgnc ID: first by ENSG and then by symbol
        my_ensg = lines[header.index("ensg_gene_name")]
        my_symbol = lines[header.index("hgnc_gene_name")]
        my_hgnc = hgnc_ids.get(my_ensg, "NA")
        if my_hgnc == "NA":
            my_hgnc = hgnc_ids.get(my_symbol, "NA")
        my_omim_id = lines[header.index("# omim_disease_id")]
        #Eliminate abbreviations from mimTitles file.
        lines[1] = lines[1].split(";")[0]
        omim_data[my_hgnc][my_omim_id] = dict(zip(header, lines))
        orphas = lines[header.index("orphanet")].split(";")
        for o in orphas:
            #Index by orphanet disease id
            orphanet_mappings[o] = my_omim_id
        #Save to extract HPO, DO terms later
        omim_data_by_disease[my_omim_id][my_hgnc] = dict(zip(header, lines))
        omim_names[my_omim_id] = lines[header.index("omim_disease_name")]

orphanet_data =  dd(lambda: dd(dict))
#keys: gene_hgnc_id, disease_omim_id (or disease name if absent)
orphanet_data_by_disease = dd(lambda: dd(dict))
with open(orphanet, 'rb') as orphanet_h:
    header = orphanet_h.readline().decode('utf8').strip().split("\t")
    for line in orphanet_h:
        lines = line.decode('utf8').strip().split("\t")
        my_hgnc = lines[header.index("hgnc_id")]
        my_orphanet_id = lines[header.index("orphanet_id")]
        my_omim_id = orphanet_mappings.get(my_orphanet_id, "NA")
        #orphanet_data
        #If no MIM disease id, use disease name
        if my_omim_id == "NA":
            disease_name = lines[header.index("disease_name")]
            orphanet_data[my_hgnc][disease_name] = dict(zip(header, lines))
        else:
            orphanet_data[my_hgnc][my_omim_id] = dict(zip(header, lines))
        #Save to extract HPO, DO terms later
        orphanet_data_by_disease[my_omim_id][my_hgnc] = dict(zip(header, lines))

decipher_data =  dd(lambda: dd(dict))
#keys: gene_hgnc_id, disease_omim_id (or disease name if absent)
decipher_data_by_disease = dd(lambda: dd(dict))
with open(decipher, 'rb') as decipher_h:
    header = decipher_h.readline().decode('utf8').strip().split("\t")
    header.append("ensg_gene_name")
    for line in decipher_h:
        lines = line.decode('utf8').strip().split("\t")
        my_hgnc = lines[header.index("hgnc_id")]
        my_omim_id = lines[header.index("omim_disease_id")]
        disease_name = lines[header.index("disease_name")]
        #Add Ensembl ID.
        my_ensg = hgnc_ids.get(my_hgnc, "NA")
        lines.append(my_ensg)
        if my_omim_id == "NA":
            decipher_data[my_hgnc][disease_name] = dict(zip(header, lines))
        else:
            decipher_data[my_hgnc][my_omim_id] = dict(zip(header, lines))
        #Save to extract HPO, DO terms later
        decipher_data_by_disease[my_omim_id][my_hgnc] = dict(zip(header, lines))

ge_data = dd(lambda: dd(dict))
#keys: gene_hgnc_id, disease_omim_id (or disease name if absent)
with open(ge, 'rb') as ge_h:
    header = ge_h.readline().decode('utf8').strip().split("\t")
    for line in ge_h:
        lines = line.decode('utf8').strip().split("\t")
        my_hgnc = lines[header.index("hgnc_id")]
        my_omim_id = lines[header.index("# omim_disease_id")]
        ge_data[my_hgnc][my_omim_id] = dict(zip(header, lines))

final_results_table = dd(lambda: dd(dict))

def add_omim(omim_data, final_results_table):
    for g in omim_data:
        for d in omim_data[g]:
            entry = omim_data[g][d]
            entry.pop('# omim_disease_id', None)
            entry.pop('cyto_location', None)
            entry.pop('orphanet', None)
            entry["source"] = {"OMIM"}
            link = "https://www.omim.org/entry/" + d
            entry["omim_link"] = link
            final_results_table[g][d] = entry
    return final_results_table  

def add_orphanet(orphanet_data, final_results_table):
    for g in orphanet_data:
        for d in orphanet_data[g]:
            entry = orphanet_data[g][d]
            entry.pop("omim_disease_id", None)
            entry.pop("disease_description", None)
            entry.pop("hgnc_id", None)
            #May want to filter on whether Assessed or not.
            entry.pop("gene_relationship_status", None)
            entry["source"] = {"Orphanet"}
            #Check for final results entry:
            if g in final_results_table:
                if d in final_results_table[g]:
                    #Add Orphanet number
                    #Possibly add HPO terms
                    final_results_table[g][d]["orphanet_id"] = entry["orphanet_id"]
                    final_results_table[g][d]["source"].add("Orphanet")
                else:
                    entry["omim_disease_name"] = entry["disease_name"]
                    entry.pop("disease_name", "")
                    final_results_table[g][d] = entry

            else:
                entry["omim_disease_name"] = entry["disease_name"]
                entry.pop("disease_name", "")
                final_results_table[g][d] = entry
    return final_results_table


def add_decipher(decipher_data, final_results_table):
    for g in decipher_data:
        for d in decipher_data[g]:
            entry = decipher_data[g][d]
            entry["disease_name"] = entry["disease_name"].lower().capitalize()
            entry.pop("organ_specificity_list", None)
            entry.pop("omim_disease_id", None)
            entry.pop("hgnc_id", None)
            entry["source"] = {"Decipher"}
            #Check for final results entry:
            if g in final_results_table:
                if d in final_results_table[g]:
                    final_results_table[g][d]["source"].add("Decipher")
                else:
                    entry["omim_disease_name"] = entry["disease_name"]
                    entry.pop("disease_name", "")
                    final_results_table[g][d] = entry

            else:
                entry["omim_disease_name"] = entry["disease_name"]
                entry.pop("disease_name", "")
                final_results_table[g][d] = entry
    return final_results_table

def add_ge(ge_data, final_results_table):
    for g in ge_data:
        for d in ge_data[g]:
            entry = ge_data[g][d]
            entry.pop("# omim_disease_id", None)
            entry.pop("hgnc_id", None)
            #Use latest genome build for Ensembl identifiers
            entry["ensg_gene_name"] = entry["ensg_gene_name_38"]
            entry.pop("ensg_gene_name_37", None)
            entry.pop("ensg_gene_name_38", None)
            entry["source"] = {"Genomics England"}
            if g in final_results_table:
                if d in final_results_table[g]:
                    final_results_table[g][d]["source"].add("Genomics England")
                else:
                    final_results_table[g][d] = entry
            else:
                final_results_table[g][d] = entry
    return final_results_table

def update_terms(hpo, do, d, data_by_disease):
    if d in data_by_disease:
        for g in data_by_disease[d]:
            if "hpo" in data_by_disease[d][g]:
                hpos = set(data_by_disease[d][g]["hpo"].split(";"))
                my_diff = hpos.difference(hpo)
                hpo.update(my_diff)
            if "do" in data_by_disease[d][g]:
                dos = set(data_by_disease[d][g]["do"].split(";"))           
                my_diff = dos.difference(do)
                do.update(my_diff)  
    return hpo, do

def add_terms(final_results_table, omim_data_by_disease, orphanet_data_by_disease, decipher_data_by_disease):
    for g in final_results_table:
        for d in final_results_table[g]:
            hpo = set(final_results_table[g][d].get("hpo", "NA").split(";"))
            do = set(final_results_table[g][d].get("do", "NA").split(";"))
            hpo, do = update_terms(hpo, do, d, omim_data_by_disease) 
            hpo, do = update_terms(hpo, do, d, orphanet_data_by_disease)
            hpo, do = update_terms(hpo, do, d, decipher_data_by_disease) 
            #If more than 1 term in hpo and do, remove "NA"
            if len(hpo) > 1:
                hpo.discard("NA")
            if len(do) > 1:
                do.discard("NA")
            hpo_string = ";".join(hpo)
            do_string = ";".join(do)
            if len(hpo_string.strip()) == 0:
                hpo_string = "NA"
            if len(do_string.strip()) == 0:
                do_string = "NA"
            final_results_table[g][d]["hpo"] = hpo_string
            final_results_table[g][d]["do"] = do_string

    return final_results_table


final_results_table = add_omim(omim_data, final_results_table)
final_results_table = add_orphanet(orphanet_data, final_results_table)
final_results_table = add_decipher(decipher_data, final_results_table)
final_results_table = add_ge(ge_data, final_results_table)
final_results_table = add_terms(final_results_table, omim_data_by_disease, orphanet_data_by_disease, decipher_data_by_disease)

out_h = open(out, 'w', encoding="utf8")
header = ["hgnc_id",  "hgnc_gene_name", "omim_gene_id", 
"ensg_gene_name", "omim_disease_id", "orphanet_id", "disease_name",
"hpo", "do", "source", "omim_link"]
out_h.write("\t".join(header) + "\n")

#Iterate over keys in dictionary of dictionaries with get, "NA" default if not found.
sorted_g = sorted(final_results_table.keys())
for hgnc_id in sorted_g:
    sorted_d = sorted(final_results_table[hgnc_id].keys())
    for omim_disease_id in sorted_d:
        current = final_results_table[hgnc_id][omim_disease_id]
        hgnc_gene_name = current.get("hgnc_gene_name", "NA")
        omim_gene_id = current.get("omim_gene_id", "NA")
        ensg_gene_name = current.get("ensg_gene_name", "NA")
        orphanet_id = current.get("orphanet_id", "NA")
        disease_name = current.get("omim_disease_name", "NA")
        omim_link = current.get("omim_link", "NA")
        if omim_disease_id.lower() == disease_name.lower():
            omim_disease_id = "NA"
        #Use OMIM disease names as default
        if omim_disease_id in omim_names:
            disease_name = omim_names[omim_disease_id]
        #Tidying up the names again
        disease_name = disease_name.replace("{", "")
        disease_name = disease_name.replace("}", "")
        disease_name = disease_name.replace("?", "")
        if disease_name.isupper():
            disease_name = disease_name.lower().capitalize()
        hpo = current.get("hpo", "NA")
        do = current.get("do", "NA")
        sources = sorted(current["source"])
        sources = ";".join(sources)
        line = [hgnc_id,  hgnc_gene_name, omim_gene_id, 
        ensg_gene_name, omim_disease_id, orphanet_id, disease_name,
        hpo, do, sources, omim_link]
        out_h.write("\t".join(line) + "\n")
out_h.close()