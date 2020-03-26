
#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, omim, orphanet, uniprot, do, infile, outfile = argv

class Generic():
    def __init__(self, source):
        self.source = source
        self.lines = []
        self.data = dd()
    def read_in(self):
        with open (self.source) as of:
            for line in of:
                lines = line.strip().split("\t")
                self.lines.append(lines)

class Omim(Generic):
    def __init__(self, source):
        super().__init__(source)

    #Index by OMIM ID
    def read_in_omim(self):
        for lines in self.lines:
            if len(lines) == 2:
                self.data[lines[0]] = lines[1]


class Orphanet(Generic):
    def __init__(self, source):
        super().__init__(source)
    #Index by Orphanet ID
    def read_in_orphanet(self):
        for lines in self.lines:
            if lines[3] != "NA":
                self.data[lines[0]] = lines[3]

class Uniprot(Generic):
    def __init__(self, source):
        super().__init__(source)
    #Index by OMIM ID
    def read_in_uniprot(self):
        for lines in self.lines:
            self.data[lines[0]] = lines[1]

class DO(Generic):
    def __init__(self, source):
        super().__init__(source)
    #Index by DO ID
    def read_in_do(self):
        for lines in self.lines:
            self.data[lines[0]] = lines[1]    

class DisTables(Generic):
    def __init__(self, source, outfile):
        super().__init__(source)
        self.outfile = outfile
    def add_descriptions(self):
        out = open(self.outfile, 'w')
        header = self.lines[0]
        out.write("\t".join(header) + "\t" + "disease_description" + "\n")
        for lines in self.lines[1:]:
            omim_disease_id = lines[4]
            orphanet_disease_id = lines[5]
            do_id = lines[8]
            omim_description = omim_o.data.get(omim_disease_id, "NA")
            orphanet_description = orphanet_o.data.get(orphanet_disease_id, "NA")
            uniprot_description = uniprot_o.data.get(omim_disease_id, "NA")
            do_description = do_o.data.get(do_id, "NA")
            my_list = [do_description, omim_description, orphanet_description, uniprot_description]
            description_to_output = max(my_list, key=len)
            out.write("\t".join(lines) + "\t" + description_to_output + "\n")
        out.close()



omim_o = Omim(omim)
omim_o.read_in()
omim_o.read_in_omim()

orphanet_o = Orphanet(orphanet)
orphanet_o.read_in()
orphanet_o.read_in_orphanet()

uniprot_o = Uniprot(uniprot)
uniprot_o.read_in()
uniprot_o.read_in_uniprot()

do_o = DO(do)
do_o.read_in()
do_o.read_in_do()

table_o = DisTables(infile, outfile)
table_o.read_in()
table_o.add_descriptions()
