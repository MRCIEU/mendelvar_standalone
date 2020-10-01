#!/usr/bin/env python
from bs4 import BeautifulSoup
from collections import defaultdict as dd
import re 
# Make our soup
with open('en_product4.xml', 'rb') as infile:
    blob = infile.read()
soup = BeautifulSoup(blob, 'lxml')

#Disease orphanumber
all_ona = dd(list)

for gna in soup.find_all(re.compile("^disorder$"), recursive=True):
	for ona in gna.find('orphacode'):
		full_tag = ona.parent
		c = full_tag.find_next_sibling('name')
		name = c.get_text()
		all_ona[ona].append(name)
		hpos = []
		dlist = c.find_next_sibling('hpodisorderassociationlist')
		for d in dlist.find_all('hpodisorderassociation'):
			my_ids = d.find("hpoid").get_text()
			my_ids = my_ids.replace("HP:", "")
			hpos.append(my_ids)
		all_ona[ona].append(hpos)

print("orphanet_id" + "\t" + "HPO_terms")
for k in all_ona:
	my_id = k
	hpos = all_ona[k][1]
	print(my_id + "\t"  + ";".join(hpos))
