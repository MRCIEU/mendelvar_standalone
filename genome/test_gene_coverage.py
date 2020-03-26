#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

##Test gene coverage of our genome interval file relative to that of our gene-disease database.
##First try to match by HGNC ID, then by ENSG ID and then by symbol.

script, my_interval, my_database = argv

interval_by_hgnc_id = dd()
interval_by_ensg_id = dd()
interval_by_symbol_id = dd()

with open(my_interval) as mh:
    for line in mh:
        lines = line.strip().split("\t")
        my_ensg = lines[6]
        my_symbol = lines[7]
        my_hgnc_id = lines[8]
        if my_ensg != "NA":
        	interval_by_ensg_id[my_ensg] = lines
        if my_hgnc_id != "NA":
        	interval_by_hgnc_id[my_hgnc_id] = lines
        if my_symbol != "NA":
        	interval_by_symbol_id[my_symbol] = lines

hgnc_success = 0
ensg_success = 0
symbol_success = 0

with open(my_database, 'rb') as my_database_h:
    header = my_database_h.readline().decode('utf8').strip().split("\t")
    for line in my_database_h:
        lines = line.decode('utf8').strip().split("\t")
        my_ensg = lines[header.index("ensg_gene_name")]
        my_symbol = lines[header.index("hgnc_gene_name")]
        my_hgnc_id = lines[header.index("hgnc_id")]
        if my_hgnc_id in interval_by_hgnc_id:
        	hgnc_success += 1
        elif my_ensg in interval_by_ensg_id:
        	ensg_success += 1
        elif my_symbol in interval_by_symbol_id:
        	symbol_success += 1
        else:
        	print(lines[0], lines[1], lines[3])
        	#print(lines)

print("We mapped %d, %d, %d gene locations by HGNC ID, ENSG ID and symbol ID, respectively." % (hgnc_success, ensg_success, symbol_success))
