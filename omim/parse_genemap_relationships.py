#!/usr/bin/env python
from GeneMap import GeneMap
from sys import argv
import re
#Print genemap relationships of the form: gene mim id, gene approved symbol, gene ensg id, condition mim number, condition name

script, genemap_f = argv

def filter_genemap (genemap):
	for line in genemap:
	#Split by condition - they are seperated by semicolon
		ids = line[12].split(";")
		for my_id in ids:
	#Save the index containing the MIM phenotype name and mapping key
			my_index = 0
	#Split by "," to access different fields for the same condition
			my_id_details = my_id.split(",")
	#Remove whitespaces
			my_id_details = [mi.strip() for mi in my_id_details]
			my_name = my_id_details[0]
			for my_i, my_id in enumerate(my_id_details):
				MyRe = r"(\d+)\s+\((\d)\)"
				MyResult = re.search(MyRe, my_id)
				if MyResult:
					my_omim_id = MyResult.group(1)
					map_key = MyResult.group(2)
					my_index = my_i
					print(line[5] + "\t" + line[8] + "\t" + line[10] + "\t" + my_omim_id + "\t"  + ", ".join(my_id_details[0:my_index]))
	

omim = GeneMap(genemap_f)		
omim.read_in_genemap()
print ("gene_mim" + "\t" + "gene_hgcn" + "\t" + "gene_ensembl" + 
	"\t" + "disease_mim" + "\t" + "disease_name")
filter_genemap(omim.lines)