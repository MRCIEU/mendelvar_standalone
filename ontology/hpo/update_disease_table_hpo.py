#!/usr/bin/env python
from sys import argv

#Update the disease-gene file generated in the disease_integration.sh scrip with additional HPO annotations present in the ontology. 
#Print any new additions!

script, hpo, disease_table, out = argv

omim_data = dict()
orpha_data = dict()
decipher_data = dict()

seen_header = 0
with open(hpo, 'r') as hpo_h:
    for line in hpo_h:
        if line.startswith("#DatabaseID"):
            headers = line.strip().split("\t")
            seen_header = 1
        elif line.startswith("#"):
            pass
        else:
            lines = line.strip().split("\t")
            #Check if not Qualifier, if so skip
            if lines[headers.index("Qualifier")] == "NOT": continue
            my_omim_id = ""
            my_orpha_id = ""
            my_decipher_name = ""
            #Check if dealing with phenotypic abnormality only
            if lines[headers.index("Aspect")] == "P":
                db = lines[headers.index("#DatabaseID")]
                hpo_term = lines[headers.index("HPO_ID")].replace("HP:", "")
                if db.startswith("OMIM:"):
                    #Extract OMIM ID
                    my_omim_id = db.replace("OMIM:", "")
                    omim_data[my_omim_id] = hpo_term
                elif db.startswith("ORPHA:"):
                    #Extract Orphanet ID
                    my_orpha_id = db.replace("ORPHA:", "")
                    orpha_data[my_orpha_id] = hpo_term
                elif db.startswith("DECIPHER:"):
                    #Extract disease title
                    my_decipher_name = lines[headers.index("DiseaseName")].lower()
                    decipher_data[my_decipher_name] = hpo_term

out_h = open(out, 'w')

with open(disease_table, 'r') as dt_h:
    header = dt_h.readline()
    headers = header.strip().split("\t")
    out_h.write(header)
    for line in dt_h:
        lines = line.strip().split("\t")
        omim_disease_id = lines[headers.index("omim_disease_id")]
        orphanet_id = lines[headers.index("orphanet_id")]
        disease_name = lines[headers.index("disease_name")].lower()
        hpos = lines[headers.index("hpo")].split(";")
        if omim_disease_id in omim_data:
            my_term = omim_data[omim_disease_id]
            if my_term not in hpos:
                    hpos.append(my_term)
                    if hpos == ["NA"]:
                        pass
        if orphanet_id in orpha_data:
            my_term2 = orpha_data[orphanet_id]
            if (my_term2 not in hpos):
                hpos.append(my_term2)            
                if hpos == ["NA"]:
                    pass

        if disease_name in decipher_data:
            my_term3 = decipher_data[disease_name]
            if (my_term3 not in hpos):
                hpos.append(my_term3)             
                if hpos == ["NA"]:
                    pass

        hpos = set(hpos)
        if (len(hpos) > 1 and "NA" in hpos):
            hpos.remove("NA")
        lines[headers.index("hpo")] = ";".join(hpos)
        out_h.write("\t".join(lines) + "\n")



out_h.close()
