#!/usr/bin/env python
from sys import argv
import re
script, do = argv

with open (do) as do_h:
	all_do = do_h.read()
	#Split by term
	terms = all_do.split("[Term]")
	for t in terms:
		temp_dict = {}
		#All the fields
		relevant_lines = set(t.split("\n"))
		for n in relevant_lines:
			if (n.startswith("id:") and "DOID:" in n):
				temp_dict["id"] = n.split("DOID:")[1]
			elif n.startswith("def:"):
				n = n.replace('"', '')
				n = n.replace("_", " ")
				n = re.sub(r"{comment.*", "", n)
				n = n.replace('def:', '')
				n = re.sub(r"\[url:.*\]", "", n)
				temp_dict["def"] = n
		if ("def" in temp_dict and "id" in temp_dict):
			print(temp_dict["id"] + "\t" + temp_dict["def"])