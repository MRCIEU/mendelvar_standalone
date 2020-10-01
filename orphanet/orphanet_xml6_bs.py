#!/usr/bin/env python
from bs4 import BeautifulSoup
import re
# Make our soup
with open('en_product6.xml', 'rb') as infile:
    blob = infile.read()
soup = BeautifulSoup(blob, 'lxml')

#Disease orphanumber
all_ona = []
#Disease name
all_dna = []
#Disease-gene association status
all_sna = []
#Gene symbol (HGCN)
all_sg = []
#Gene ENSG ID
all_ensg = []
#Gene HGNC ID
all_hna = []
#Gene omim ID
all_omim = []

def find_identifier(d, database):
    s = d.find_all("source", string=re.compile(database))
    if s:
        a = s[0].find_next_sibling('reference')
        out = (a.get_text())
    else:
        out = "NA"
    return out

#Get Orphanumber 
for gna in soup.find_all(re.compile("^disorder$"), recursive=True):
    for ona in gna.find('orphacode'):
        full_tag = ona.parent
        c = full_tag.find_next_sibling('name')
        name = c.get_text()
        dlist = c.find_next_sibling('disordergeneassociationlist')
        for d in dlist.find_all('disordergeneassociation'):
            all_ona.append(ona)
            all_dna.append(name)
            #Symbol
            all_sg.append(d.find("symbol").get_text())
            #ENSG
            all_ensg.append(find_identifier(d, "Ensembl"))
            #HGNC
            all_hna.append(find_identifier(d, "HGNC"))
            #OMIM
            all_omim.append(find_identifier(d, "OMIM"))
            status = d.find('disordergeneassociationstatus')
            status_name = status.find('name')
            all_sna.append((status_name.get_text()))
    
#zip everything together
print("orphanet_id" + "\t" + "gene_symbol" + "\t" + "gene_ensg" +
"\t" + "gene_hgnc_id" + "\t" + "gene_omim_id" + "\t" + "gene_relationship_status")
zipped = zip(all_ona, all_sg, all_ensg, all_hna, all_omim, all_sna)
for my_t in zipped:
    print("\t".join(my_t))
