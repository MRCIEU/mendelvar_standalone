#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, my_input, my_31_input, hgnc = argv

hgnc_ids = dd()
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

ver_31_mappings = dd()
#Match hgnc symbols to identifiers
with open(my_31_input) as mh31:
    for line in mh31:
        lines = line.strip().split("\t")
        my_ensg = lines[6]
        my_symbol = lines[7]
        my_hgnc_id = lines[8]
        ver_31_mappings[my_symbol] = my_hgnc_id
        #Check if hgnc ID matches symbol.


def get_hgnc(my_ensg, my_symbol):
    if my_ensg in hgnc_ids:
        my_hgnc_id = hgnc_ids[my_ensg]
        if my_hgnc_id in hgnc_ids:
            if my_symbol == hgnc_ids[my_hgnc_id]:
                pass
            else:
                my_symbol = hgnc_ids[my_hgnc_id]
        else:
            pass
    else:
        if my_symbol in hgnc_ids:
            my_hgnc_id = hgnc_ids[my_symbol]
        else:
            my_hgnc_id = "NA"
    return my_symbol, my_hgnc_id

with open(my_input) as mh:
    for line in mh:
        lines = line.strip().split("\t")
        my_ensg = lines[6]
        my_symbol = lines[7]
        my_hgnc_id = lines[8]
        lines[7], lines[8] = get_hgnc(my_ensg, my_symbol)
        #Check if my symbol present in mappings to hgnc identifiers from v31.
        """
        if my_symbol in ver_31_mappings:
            my_hgnc_id = ver_31_mappings[my_symbol]
            #Swap to HGNC id.
            lines[8] = my_hgnc_id 
            if my_hgnc_id in hgnc_ids:
                #Check if hgnc ID matches symbol.
                if my_symbol == hgnc_ids[my_hgnc_id]:
                    pass
                else:
                    #Swap to the name in the HUGO file.
                    lines[7] = hgnc_ids[my_hgnc_id]
            else:
                #Check if we can get hgnc id and symbol through ENSG 
                lines[7], lines[8] = get_hgnc(my_ensg, my_symbol)
                
        else:
            #Check if we can get hgnc id and symbol through ENSG 
            lines[7], lines[8] = get_hgnc(my_ensg, my_symbol)
        """
        print("\t".join(lines))
            