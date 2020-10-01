#!/usr/bin/env python
from bs4 import BeautifulSoup
# Make our soup
with open('en_product1.xml', 'rb') as infile:
    blob = infile.read()
# Use LXML for blazing speed
soup = BeautifulSoup(blob, 'lxml')

all_ona = []
all_dna = []
all_cna = []
all_omim = []

#Get Orphanumber
for ona in soup.find_all('orphacode'):
	my_p = ona.parent.find_all('expertlink')
	if my_p:
		all_ona.append(ona.get_text())

#Get disease name
for dna in soup.find_all('name'):
	my_p = dna.parent.find_all('expertlink')
	if my_p:
		all_dna.append(dna.get_text())

#Get disease description
for d in soup.find_all('disorder'):
	if d.find('contents'):
		all_cna.append(d.find('contents').get_text())
	else:
		all_cna.append("NA")

#Get OMIM ID
for d in soup.find_all('disorder'):
	if d.find_all('source'):
		source = d.find_all("source")
		omim_flag = 0
		for s in source:
			if s.get_text().strip() == "OMIM":
				all_omim.append(s.parent.find("reference").get_text())
				omim_flag = 1
				break
		if not omim_flag:
			all_omim.append("NA")
	else:
		all_omim.append("NA")
			
#zip everything together
print("orphanet_id" + "\t" + "disease_omim_id" + "\t" + "disease_name" + "\t" + "disease_description")
zipped = zip(all_ona, all_omim, all_dna, all_cna)
for my_t in zipped:
	print("\t".join(my_t))
