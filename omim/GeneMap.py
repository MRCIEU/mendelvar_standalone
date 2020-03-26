#!/usr/bin/env python
import re

class GeneMap:
	def __init__(self, genemap_f):
		self.genemap = genemap_f
		self.lines = []
		self.final_parsed = []
		self.omim_disease = set()

	def read_in_genemap(self):
		with open(self.genemap) as gf:
			for line in gf:
				if line.startswith("#"):
					next
				else:
					lines = line.strip().split("\t")
					#Check if enough fields to contain phenotype and gene name.
					#Field 10 checks for ENSG identifier and field 12 for phenotype.
					if (len(lines) > 12 and lines[10]!="" and lines[12]!=""):
						self.lines.append(lines)
			return self.lines

				
	def filter_genemap(self):
		for line in self.lines:
			#Split by condition - they are seperated by semicolon
			ids = line[12].split(";")
			for my_id in ids:
				#Save the index containing the MIM phenotype name and mapping key
				my_index = 0
				#Split by "," to access different fields for the same condition
				my_id_details = my_id.split(",")
				#Remove whitespaces
				my_id_details = [mi.strip() for mi in my_id_details]
				my_name = my_id_details[0]
				#Eliminate non-diseases in []
				if my_name.startswith(r'['):
					next
				else:
					#Check which field contains condition ID
					for my_i, my_id in enumerate(my_id_details):
						MyRe = r"(\d+)\s+\((3)\)"
						MyResult = re.search(MyRe, my_id)
						if MyResult:
							my_omim_id = MyResult.group(1)
							map_key = MyResult.group(2)
							my_index = my_i
							#Print result
							#Omim_disease_id, omim_disease_name, omim_gene_id, 
							#hgnc_name, ensg_name, cyto_location	
							self.final_parsed.append([my_omim_id, ", ".join(my_id_details[0:my_index]),  
							line[5], line[8], line[10], line[3]])
							self.omim_disease.add(my_omim_id)

		return self.final_parsed

	#Same as above but any mapping key possible
	def filter_genemap_no_mk(self):
		for line in self.lines:
			#Split by condition - they are seperated by semicolon
			ids = line[12].split(";")
			for my_id in ids:
				#Save the index containing the MIM phenotype name and mapping key
				my_index = 0
				#Split by "," to access different fields for the same condition
				my_id_details = my_id.split(",")
				#Remove whitespaces
				my_id_details = [mi.strip() for mi in my_id_details]
				my_name = my_id_details[0]
				#Eliminate non-diseases in []
				if my_name.startswith(r'['):
					next
				else:
					#Check which field contains condition ID
					for my_i, my_id in enumerate(my_id_details):
						MyRe = r"(\d+)\s+\((\d)\)"
						MyResult = re.search(MyRe, my_id)
						if MyResult:
							my_omim_id = MyResult.group(1)
							map_key = MyResult.group(2)
							my_index = my_i
							#Print result
							#Omim_disease_id, omim_disease_name, omim_gene_id, 
							#hgnc_name, ensg_name, cyto_location	
							self.final_parsed.append([my_omim_id, ", ".join(my_id_details[0:my_index]),  
							line[5], line[8], line[10], line[3]])
							self.omim_disease.add(my_omim_id)

		return self.final_parsed

	def print_genemap(self):
		print("# omim_disease_id" + "\t" +"omim_disease_name" + "\t"
		+ "omim_gene_id" + "\t" + "hgnc_gene_name" + "\t" + "ensg_gene_name" + "\t" + "cyto_location")
		for my_list in self.final_parsed:
			print("\t".join(my_list))
		return 1

	def update_final_parsed(self, to_keep):
		temp_list = []
		for my_index, my_list in enumerate(self.final_parsed):
			if my_list[0] in to_keep:
				temp_list.append(my_list)
				
		self.final_parsed = temp_list
		return self.final_parsed


