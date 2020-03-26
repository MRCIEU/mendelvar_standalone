#!/usr/bin/env python
from sys import argv
import re
#Creates two outputs based on mimTitles.txt: 
#-File with all disease IDs currently in OMIM (as marked by Plus, Percent, Null, Number Sign)
##-File with all gene IDs currently in OMIM (as marked by Plus, Asterisk)
#-file which maps deprecated mim ids to the current ones (as marked by Caret)
script, mimtitles = argv

caret = open('omim_moved.txt', 'w')
all_disease = open('mimTitles_processed.txt' ,'w')
all_genes = open('mimTitles_genes.txt' ,'w')


with open (mimtitles) as mimt:
	for line in mimt:
		if not line.startswith("#"):
			lines = line.strip().split("\t")
			if lines[0].strip() == "Caret":
				old_id = lines[1]
				new_id = lines[2].replace("MOVED TO ", "").strip()
				#Write out old id followed by new id
				if new_id == "REMOVED FROM DATABASE":
					pass
				else:
					new_id = new_id.split()[0]
					new_id = re.sub(r"\D*(\d+)\D*", r"\1", new_id)
				caret.write(old_id + "\t" + new_id + "\n")
				#In certain cases, we have "REMOVED FROM DATABASE", ID enclosed in curly braces and two ids joined by AND.
				#In the first case, include in the output (so that e.g. Decipher entries with those MIM ID are known)
				#In the latter case, pick up only the first ID (split on space).
				#In the second case, get rid of curly braces (match only digits)
			if (lines[0].strip() != "Asterisk" and lines[0].strip() != "Caret"):
				all_disease.write(lines[1] + "\n")
			if (lines[0].strip() == "Asterisk" or lines[0].strip() == "Plus"):
				all_genes.write(lines[1] + "\n")


caret.close()
all_disease.close()
all_genes.close()