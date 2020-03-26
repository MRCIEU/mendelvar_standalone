#!/usr/bin/env python
import glob, json
from collections import defaultdict as dd

output_dict = dd(list)
second_dict = dd(list)
for filename in glob.iglob('json/*.json'):
	with open (filename, 'r') as db_h:
		my_id = None
		response_dict = json.load(db_h)
		disease_id = response_dict["omim"]["entryList"][0]["entry"]["mimNumber"]
		if ("phenotypicSeriesExists" in response_dict["omim"]["entryList"][0]["entry"] and response_dict["omim"]["entryList"][0]["entry"]["phenotypicSeriesExists"] == True):
			if ("phenotypeMapExists" in response_dict["omim"]["entryList"][0]["entry"] and response_dict["omim"]["entryList"][0]["entry"]["phenotypeMapExists"] == True):
				my_id = response_dict["omim"]["entryList"][0]["entry"]["phenotypeMapList"][0]["phenotypeMap"].get("phenotypicSeriesNumber", None)
			elif ("geneMapExists" in response_dict["omim"]["entryList"][0]["entry"] and response_dict["omim"]["entryList"][0]["entry"]["geneMapExists"] == True):
				my_id = response_dict["omim"]["entryList"][0]["entry"]["geneMap"]["phenotypeMapList"][0]["phenotypeMap"].get("phenotypicSeriesNumber", None)
		if my_id:
			my_ids = my_id.strip().split(",")
			my_ids = [i.replace("PS", "") for i in my_ids]
			[output_dict[i].append(disease_id) for i in my_ids]

for v in output_dict.values():
	print(v)
	for my_i, my_e in enumerate(v):
		second_dict[my_e].extend(v[0:my_i])
		second_dict[my_e].extend(v[my_i+1:])
#Output dictionary with keys: phenotype series IDs, values: disease MIMs belonging to the given Phenotype Series.
with open('phenotype_series.json', 'w') as f:
	json.dump(output_dict, f, indent=4, sort_keys=True)

with open('phenotype_series2.json', 'w') as f2:
	json.dump(second_dict, f2, indent=4, sort_keys=True)

