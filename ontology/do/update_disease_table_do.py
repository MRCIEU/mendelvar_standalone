#!/usr/bin/env python
from sys import argv

#Update the disease-gene file generated in the disease_integration.sh scrip with additional DO annotations present in the ontology. 
#Print any new additions!

script, do, disease_table, out = argv

omim_data = dict()
orpha_data = dict()

with open(do, 'r') as do_h:
    all_do = do_h.read()
    terms = all_do.split("[Term]")
    terms = [t.strip() for t in terms]
    for t in terms:
        #Reset id
        my_id = ""
        my_omim = ""
        my_orpha = ""
        for line in t.split("\n"):
            if line.startswith("id:"):
                my_id = line.replace("id: DOID:", "")
            elif line.startswith("xref: OMIM:"):
                my_omim = line.replace("xref: OMIM:", "")
            elif line.startswith("xref: ORDO:"):
                my_orpha = line.replace("xref: ORDO:", "")
        if (my_id and my_omim):
            omim_data[my_omim] = my_id
        if (my_id and my_orpha):
            orpha_data[my_orpha] = my_id

out_h = open(out, 'w')

with open(disease_table, 'r') as dt_h:
    header = dt_h.readline()
    headers = header.strip().split("\t")
    out_h.write(header)
    for line in dt_h:
        lines = line.strip().split("\t")
        omim_disease_id = lines[headers.index("omim_disease_id")]
        orphanet_id = lines[headers.index("orphanet_id")]
        dos = lines[headers.index("do")].split(";")
        if omim_disease_id in omim_data:
            my_term = omim_data[omim_disease_id]
            if my_term not in dos:
                dos.append(my_term)
                if dos == ["NA"]:
                    pass
        if orphanet_id in orpha_data:
            my_term2 = orpha_data[orphanet_id]
            if (my_term2 not in dos):
                dos.append(my_term2)            
                if dos == ["NA"]:
                    pass
        #Get only unique DO terms
        dos = set(dos)
        if (len(dos) > 1 and "NA" in dos):
            dos.remove("NA")
        lines[headers.index("do")] = ";".join(dos)
        out_h.write("\t".join(lines) + "\n")


out_h.close()