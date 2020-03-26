#!/usr/bin/env python
from sys import argv
import json, glob, re
script, mimtitles_p = argv
#Input: file with all disease IDS in OMIM:
#Read into a set
mim = set()

with open(mimtitles_p) as mph:
	for line in mph:
		mim.add(int(line.strip()))

#Read in all downloaded JSON files.
for filename in glob.iglob('json/*.json'):
	with open (filename, 'r') as db_h:
		response_dict = json.load(db_h)
		disease_id = response_dict["omim"]["entryList"][0]["entry"]["mimNumber"]
		if disease_id in mim:
			text = response_dict["omim"]["entryList"][0]["entry"].get("textSectionList")
			if text:
				for t in response_dict["omim"]["entryList"][0]["entry"]["textSectionList"]:
					text2 = t["textSection"]["textSectionName"]
					if text2 == "description":
						description = t["textSection"]["textSectionContent"]
						#Cut down everything after Subheadings
						res = description.partition("<Subhead>")
						if len(res[0]) > 2:
							res = res[0]
						else:
							res = res[1]
						#Cut down everything after "For" review etc. links
						res2 = res.partition("For")[0].split("\n")[0].replace("{", "").replace("}", "")
						#Remove additional references 
						res2 = re.sub(r"\d+:", "", res2)
						#Print disease omim ID and description
						print(str(disease_id) + "\t" + res2)