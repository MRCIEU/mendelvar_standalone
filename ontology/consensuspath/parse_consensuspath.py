#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, pathdb, hgnc, out = argv

hgnc_ids = dict()
with open(hgnc, 'rb') as hgnc_h:
    header = hgnc_h.readline().decode('utf8').strip().split("\t")
    for line in hgnc_h:
        lines = line.decode('utf8').strip().split("\t")
        hg_id = lines[0].replace("HGNC:", "")
        hg_symbol = lines[1]
        hgnc_ids[hg_symbol] = hg_id

out_h = open(out, "w")

with open(pathdb) as pathdb_h:
    header = pathdb_h.readline()
    counter = 0
    for line in pathdb_h:
        lines = line.strip().split("\t")
        title = lines[0]
        pathway = lines[1]
        if pathway == "None":
            counter += 1 
            pathway = "custom_id_" + str(counter)
        symbols = lines[3].strip().split(",")
        for s in symbols:
            my_id = hgnc_ids.get(s, "NA")
            if s == "NA":
                print ("Error %s" % s)
                continue
            out_h.write(my_id + "\t" + pathway + "\t" + title + "\n")


out_h.close()