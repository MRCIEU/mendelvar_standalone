#!/usr/bin/env python
from sys import argv
#Create a file to be supplied to the omim_api.py file which includes all the OMIM diseases with
#no known molecular basis (i.e. Plus, Percent, Null). This is so that we can extract HPO, DO
#for these diseases which could be then annotated by genes from another database.
script, mimtitles = argv

mim_out = open('mimtitles_parsed_master.txt', 'w')

mim_out.write("# omim_disease_id" + "\t" +"omim_disease_name" + "\t"
+ "omim_gene_id" + "\t" + "hgnc_gene_name" + "\t" + "ensg_gene_name" + "\t" + "cyto_location" + "\n")

with open (mimtitles) as mimt:
	for line in mimt:
		if not line.startswith("#"):
			lines = line.strip().split("\t")
			if (lines[0].strip() != "Asterisk" and lines[0].strip() != "Caret" and lines[0].strip() != "Number Sign"):
				mim_out.write(lines[1] + "\t" + lines[2] + "\t" + "NA" +
				"\t" + "NA" + "\t" + "NA" + "\t" + "NA" + "\n")


mim_out.close()