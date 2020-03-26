#!/usr/bin/env python
from sys import argv

script, disease_table, out = argv

out_h = open(out, 'w')
out_h.write("!gaf-version: 2.1" + "\n")

with open(disease_table, 'r') as dt_h:
    header = dt_h.readline()
    headers = header.strip().split("\t")
    db = "HGNC"    
    qualifier = ""
    for line in dt_h:
        lines = line.strip().split("\t")
        db_object_id = lines[headers.index("hgnc_id")]
        db_symbol = lines[headers.index("hgnc_gene_name")]
        go_ids = lines[headers.index("hpo")]     
        db_ref = lines[headers.index("source")]
        evidence_code = "TAS"
        with_from = "NA"
        aspect = "M"
        db_object_name = ""
        db_object_synonym = ""
        db_object_type = "protein"
        taxon = "taxon:9606"
        date = "20090118"
        assigned = "MKS"
        ext = ""
        gene_product = ""
        if go_ids != "NA":
            all_go_ids = go_ids.split(";")
            for a in all_go_ids:
                my_line = [db, db_object_id, db_symbol, qualifier, a, 
                db_ref, evidence_code, with_from, aspect, db_object_name,
                db_object_synonym, db_object_type, taxon, date, assigned,
                ext, gene_product]
                out_h.write("\t".join(my_line) + "\n")


out_h.close()