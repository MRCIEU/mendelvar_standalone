#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, my_input, hgnc = argv

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
            hgnc_ids[hg_id] = hg_symbol
            hgnc_ids[hg_symbol] = hg_id
        else:
            hgnc_ids[hg_symbol] = hg_id
            hgnc_ids[hg_id] = hg_symbol


with open(my_input) as mh:
    for line in mh:
        lines = line.strip().split("\t")
        my_ensg = lines[6]
        my_symbol = lines[7]
        my_hgnc_id = lines[8]
        #Check if hgnc ID matches symbol.
        if my_hgnc_id in hgnc_ids:
            if my_symbol == hgnc_ids[my_hgnc_id]:
                pass
            else:
                #Swap to the name in the HUGO file.
                lines[7] = hgnc_ids[my_hgnc_id]
        else:
            pass

        print("\t".join(lines))
