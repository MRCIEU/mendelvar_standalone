#!/usr/bin/env python
import argparse
import statsmodels.stats.multitest as smm
'''
INRICH WRAPPER
parses INRICH output for display
adds HUGO gene names and ENSG identifiers matching a given gene dataset
add FDR values
#
'''

class Inrich:
	def __init__(self, inrich):
		self.inrich = inrich
		self.my_data = dict()
		self.my_pvalues = list()
	def read_in_inrich(self):
		"""
		Read in data from INRICH output. Read into a dictionary hashed by GO term.
		"""
		my_marker = 0
		try:
			with open(self.inrich, 'r') as inrich_fh:
				for line in inrich_fh:
					if (line.startswith("_") or not line.strip()):
						my_marker = 0
					if line.startswith("T_Size"):
						line = next(inrich_fh)
						line = next(inrich_fh)
						my_marker=1
					if my_marker == 1:
						lines = line.strip().split()
						self.my_data[lines[4]] = lines[0:4]
						self.my_pvalues.append(float(lines[2]))
				return self
		except FileNotFoundError:
			print ("Warning: Inrich file not found")

	def pvalue(self):
		return self.add_pvalues(self.my_pvalues)

	@staticmethod
	def add_pvalues(p):
		"""
		Create a dictionary with empirical p-values as keys and FDR-corrected p-values as values.
		"""
		fdr_pval = list(smm.multipletests(p, method='fdr_bh', is_sorted=False, returnsorted=False)[1:2][0])
		zipped_dict = dict(zip(p, fdr_pval))
		return zipped_dict


class GeneMap:
	def __init__(self, gene_map):
		self.gene_map = gene_map
	def read_in_gene_map(self):
		"""
		Read in gene map. Read into a dictionary hashed by gene id.
		"""
		self.my_data = dict()
		try:
			with open(self.gene_map, 'r') as gene_map_fh:
				for line in gene_map_fh:
					lines = line.strip().split()
					self.my_data[lines[3]] = lines[4]
				return self.my_data
		except FileNotFoundError:
			print ("Warning: GeneMap file not found")

class MatchingGenes:
	def __init__(self, matching_genes):
		self.matching_genes = matching_genes
	def read_in_genes(self):
		"""
		Read in genes overlapping the interval. Read into a dictionary hashed by gene id, with values being HGNC ID gene names.
		"""
		self.my_data = dict()
		try:
			with open(self.matching_genes, 'r') as matching_genes_fh:
				for line in matching_genes_fh:
					lines = line.strip().split()
					self.my_data[lines[5]] = lines[4]
				return self.my_data
		except FileNotFoundError:
			print ("Warning: MatchingGenes file not found")

class GeneSet:
	def __init__(self, gene_set):
		self.gene_set = gene_set
	def read_in_gene_set(self):
		"""
		Read in gene set. Read into a dictionary hashed by gene set id.
		"""
		self.my_data = dict()
		try:
			with open(self.gene_set, 'r') as gene_set_fh:
				for line in gene_set_fh:
					lines = line.strip().split("\t")
					if lines[1] in self.my_data:
						gene_list = self.my_data[lines[1]][1]
						gene_list.add(lines[0])
						self.my_data[lines[1]] = [lines[2], gene_list]
					else:
						gene_list = set()
						gene_list.add(lines[0])
						self.my_data[lines[1]] = [lines[2], gene_list]
				return self.my_data
		except FileNotFoundError:
			print ("Warning: GeneSet file not found")

class PrintOutput:
	def __init__(self, inrich_in,inrich_pval, gene_map, gene_set, out):
		self.inrich_in = inrich_in.my_data
		self.inrich_pval = inrich_pval
		self.gene_map = gene_map
		self.gene_set = gene_set
		self.out = out

	def combine_output(self):
		with open(self.out, "w") as output:
			output.write("Gene_set_size" + "\t" + "Overlap_size" + "\t" + "Overlapping_genes"
			+ "\t" + "Empirical_p_value"+ "\t" + "Bootstrap_p_value" + "\t" + "FDR_p_value" 
			+ "\t" + "Gene_set_ID" + "\t" + "Gene_set_description" + "\n")
			#Sort INRICH output by empirical p-value
			my_sorted = sorted(self.inrich_in.items(), key=lambda x: x[1][2])
			for key, values in my_sorted:
				output.write(values[0] + "\t" + values[1] + "\t" + "NA"
				+ "\t" + values[2] + "\t" + values[3] + "\t" + str(self.inrich_pval[round(float(values[2]),6)])
				+ "\t" + key + "\t" + self.gene_set[key][0] + "\n")
