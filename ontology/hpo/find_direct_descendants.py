#!/usr/bin/env python
from sys import argv
import re

script, full_inrich, full_obo, to_keep, out = argv

#Read a list of IDs to find direct descendants of
root = set()
with open(to_keep) as ks:
    for line in ks:
        root.add(line.strip())

descendants = set()
with open(full_obo, 'r') as full_obo_h:
    all_do = full_obo_h.read()
    terms = all_do.split("[Term]")
    terms = [t.strip() for t in terms]
    for t in terms:
        my_id=""
        for line in t.split("\n"):
            if line.startswith("id:"):
                match = re.search(r'\d+', line)
                my_id = match.group()
            if line.startswith("is_a:"):
                match = re.search(r'\d+', line)
                my_des = match.group()
                #Find ID which the terms is descendent from.
                if my_des in root:
                    descendants.add(my_id)

out_fh = open(out, 'w')
print(descendants)
with open(full_inrich, 'r') as full_inrich_h:
    for line in full_inrich_h:
        lines = line.strip().split()
        if lines[1] in descendants:
            out_fh.write(line)

out_fh.close()