#!/usr/bin/env python
from sys import argv
import re
from GeneMap import GeneMap

#Script to check which entries need to be rerun with our omim_api.py script to be added to the database.
script, genemap_f, omim_update = argv

class UpdateHTML:
	def __init__(self, omim_update):
		self.htmlfile = omim_update
		self.potential_targets = set()

	def filter_html(self, all_omim_ids):
		id_regex = r'/entry/(\d+)">'
		with open(self.htmlfile) as hf:
			for line in hf:
				MyResult = re.search(id_regex, line)
				if MyResult:
					 my_id = MyResult.group(1)
					 self.potential_targets.add(my_id)

		return self.potential_targets

#Basic filtering applied as normal to GeneMap files
omim = GeneMap(genemap_f)
omim.read_in_genemap()
omim.filter_genemap()
all_omim_ids = omim.omim_disease
#Read in HTML files and return potential candidates 
#for searching for disease ids.
my_html = UpdateHTML(omim_update)
to_keep = my_html.filter_html(all_omim_ids)

#Only output those lines of parsed GeneMap file whose 
#disease ids were found in the update HTML file.
omim.update_final_parsed(to_keep)
omim.print_genemap()