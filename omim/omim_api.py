#!/usr/bin/env python

import requests, json, re
from sys import argv
from collections import defaultdict as dd
from os.path import exists

#Omim output from parse_genemap.py script.
omim_in=argv[1]
#Are we carrying out an update? In that case, re-download JSON files even if they exist.
to_update=int(argv[2])

#User defined prefix for output
file_prefix = argv[3]

#Switch on if want to disable regex check 
override_regex = int(argv[4])

#Get all attributes excluding references and name of contributors
#Change the value of apiKey in search_string to your own API key
search_string = 'https://api.omim.org/api/entry?mimNumber=mkesobczyk&include=all&exclude=reference,referenceList,contributors,creationDate,editHistory&format=json&apiKey=my_api_key'

#Count number of disease-gene associations excluded due to coming from GWAS
gwas_counter = 0
#Dump the list of MIM entries that contain the GWAS keywords
gwas_list = set()
#Track MIM numbers and their associated unsuccessful http request codes
error_list = []

def _flatten_items(items, sep, prefix):
  _items = []
  for key, value in items:
    _prefix = "{}{}".format(prefix, key)
    if isinstance(value, list):
      _items.extend(_flatten_items(list(enumerate(value)), sep=sep,
                    prefix=_prefix+sep))
    elif isinstance(value, dict):
      _items.extend(_flatten_items(value.items(), sep=sep,
                    prefix=_prefix+sep))
    else:
      _items.append((_prefix, value))
  return _items


def flatten_dict(d, sep='_'):
  return dict(_flatten_items(d.items(), sep=sep, prefix=""))


def test_gwas(response_dict):
	mygwas_re = r"genome.?wide association"
	mygwas_re2 = r"gwas"
	mycase_re = r"case.?control"
	flattened = flatten_dict(response_dict)
	for my_v in flattened.values():
		#Match the 3 different regex
		MyResult = re.search(mygwas_re, str(my_v).lower())
		MyResult2 = re.search(mygwas_re2, str(my_v).lower())
		MyResult3 = re.search(mycase_re, str(my_v).lower())
		if (MyResult or MyResult2 or MyResult3):
			global gwas_counter
			gwas_counter += 1
			gwas_list.add(flattened["omim_entryList_0_entry_mimNumber"])
			print ("Number of entries eliminated due to GWAS keywords: %d" % gwas_counter)
			return 0
		else:
			pass
	return 1


#HPO regex
hpo_re = r"HP\:\d+"

#Container for the initial list of potential OMIM hits to be included in the database.
gwas = open (file_prefix + "_gwas_to_include_omim.txt", 'w')
gwas.write("# omim_disease_id" + "\t" +"omim_disease_name" + "\t"
+ "omim_gene_id" + "\t" + "hgnc_gene_name" + "\t" + "ensg_gene_name" + "\t" + "cyto_location" "\t" + "hpo" + "\t" + "do" + "\t" + "orphanet" + "\n")

dcomplex = 0
with open (omim_in) as oh:
	for line in oh:
		if not line.startswith("#"):
			lines = line.strip().split("\t")
			#Check if susceptibility phenotype
			if "{" in lines[1]:
				dcomplex = 1		
			else:
				dcomplex = 0
			#Include in the final database?
			to_include=0
			my_id=lines[0]
			#Check if we already have JSON on file.
			if (not exists("./json/" + my_id +".json") or to_update):
				print ("Downloading OMIM id: %s" %my_id)
				search_string_rep = search_string.replace("mkesobczyk", my_id)
				response = requests.get(search_string_rep)
			#Check if correct response status:
				if response.status_code != 200:
				#Save those for later to inspect
				#raise ConnectionError("Oops... something went wrong, status code is: %s" % response.status_code)
				#Print a list of omim ids and their error codes
					error_list.append((my_id,response.status_code))
					print("oops")
					continue
				else:
					response_dict = json.loads(response.content)
				#Dump raw data to file
					out = open("./json/" + my_id+".json", 'w')
					file_dump = json.dump(response_dict, out, indent=4, sort_keys=False)
					out.close()
			with open ("./json/" + my_id+".json", 'r') as db_h:
				response_dict = json.load(db_h)
			#Name of the disease
			disease_name = response_dict["omim"]["entryList"][0]["entry"]["titles"]["preferredTitle"]
			if override_regex:
				to_include = 1
			elif "FAMILIAL" in disease_name.upper():
				to_include = 1
			elif "SUSCEPTIBILITY" in disease_name.upper():
				to_include = test_gwas(response_dict)
			elif dcomplex == 1:
				to_include = test_gwas(response_dict)
			else:
				to_include = 1
			#If included, further scan for HPO and DO terms:
			if to_include:
				gwas.write("\t".join(lines))
				clinical_synposis_status = response_dict["omim"]["entryList"][0]["entry"]["clinicalSynopsisExists"]
				#Search for HPO terms
				if clinical_synposis_status:
					clinical_synopsis = response_dict["omim"]["entryList"][0]["entry"]["clinicalSynopsis"]
					all_matches_hpo = list()
					for v in clinical_synopsis.values():
						all_matches_hpo.extend(re.findall(hpo_re, str(v)))
					#Unique HPO terms
					all_matches_hpo = [a.replace("HP:", "") for a in all_matches_hpo]
					all_matches_hpo = set(all_matches_hpo)	
					gwas.write("\t" + ";".join(all_matches_hpo))
				else:
					gwas.write("\t" + "NA")
				diseaseontology_ids = response_dict["omim"]["entryList"][0]["entry"]["externalLinks"].get("diseaseOntologyIDs", "NA")
					#Check if multiple DO terms in a list
				if isinstance(diseaseontology_ids, list):
					diseaseontology_ids = ";".join(diseaseontology_ids)
					gwas.write("\t" + diseaseontology_ids)
				else:
					gwas.write("\t" + diseaseontology_ids)
				#orphanet_disease_IDs.
				orphanet_ids = response_dict["omim"]["entryList"][0]["entry"]["externalLinks"].get("orphanetDiseases", "NA")
				my_list = orphanet_ids.strip().split(";;;")
				os = [o.strip().split(";;")[0] for o in my_list]
				gwas.write("\t" +  ";".join(os))
				gwas.write("\n")

gwas.close()

filt = open (file_prefix + "_gwas_filtered_omim.txt", 'w')
[filt.write(str(x) + "\n") for x in gwas_list]
filt.close()

error = open (file_prefix + "_error_entries.txt", 'w')
[error.write(str(x[0]) + "\t" + str(x[1]) + "\n") for x in error_list]
error.close()