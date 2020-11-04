#!/usr/bin/env python
import argparse, hashlib, gzip, re

#Parse Clinvar tabbed dataset. Output two files for hg19 and hg 38.
#Fields needed: dbSNP_dbvar_id (optional), cyto_location, chromosome, start, end, allele1, allele2, effect (pathogenic etc.), hgsv notation, gene (not obligatory), disease name, disease ID - omim and orphanet, database source, database link, VCV, RCV, quality rating (stars).

#['#AlleleID', 'Type', 'Name', 'GeneID', 'GeneSymbol', 'HGNC_ID', 'ClinicalSignificance', 'ClinSigSimple', 'LastEvaluated', 'RS# (dbSNP)', 'nsv/esv (dbVar)', 'RCVaccession', 'PhenotypeIDS', 'PhenotypeList', 'Origin', 'OriginSimple', 'Assembly', 'ChromosomeAccession', 'Chromosome', 'Start', 'Stop', 'ReferenceAllele', 'AlternateAllele', 'Cytogenetic', 'ReviewStatus', 'NumberSubmitters', 'Guidelines', 'TestedInGTR', 'OtherIDs', 'SubmitterCategories', 'VariationID']

ap = argparse.ArgumentParser()
ap.add_argument('--clinvar',required=True,type=str,help='Original gzipped tabbed ClinVar file')
ap.add_argument('--md5',required=True,type=str,help='File with md5 checksum')
ap.add_argument('--omim',required=True,type=str,help='OMIM MorbidMap')
ap.add_argument('--uncertain',required=True,type=int,help='Do we filter out uncertain sig variants? Ie. keep only pathogenic, likely, pathogenic, or risk variants. 0/1 flag')
ap.add_argument('--gene',required=True,type=int,help='Do we filter out variants spanning multiple genes? 0/1 flag')
args = ap.parse_args()

class Clinvar:
	empty_field="NA"
	omim_re = r"OMIM\:(\d+)"
	orphanet_re = r"Orphanet\:(?:ORPHA)?(\d+)"
	omim_re_mm = r"(\d+) \(\d\)"

	def __init__(self, args):
		self.clinvar = args.clinvar
		self.uncertain = args.uncertain
		self.gene = args.gene
		self.md5 = args.md5
		self.omim = args.omim
		self.omim_dict = {}
		self.added_names = set()
	def check_if_corrupted(self):
		#Read in MD5 checksum
		correct_md5 = open(self.md5, "r").readline().split()[0]
		clinvar_input = open(self.clinvar, 'rb')
		current_md5 = hashlib.md5(clinvar_input.read()).hexdigest()
		assert current_md5 == correct_md5, ("Oops, our ClinVar download was corrupted!")

	def initialise_output(self):
		self.hg19 = open("hg19_clinvar.txt", 'w')
		new_header = ["dbsnp_dbvar_id", "cyto_location", "chromosome", "start", "end",
		"ref_allele", "alt_allele", "effect", "HGSV_notation", "gene", "disease_name(s)", "disease_omim", "disease_orphanet",
		"database_source", "database_link", "VCV", "RCV", "allele_id", "quality_rating"]
		self.hg19.write("\t".join(new_header) + "\n")
		self.hg19_bed = open("hg19_clinvar.bed", 'w')
		self.hg38 = open("hg38_clinvar.txt", 'w')
		self.hg38.write("\t".join(new_header) + "\n")
		self.hg38_bed = open("hg38_clinvar.bed", 'w')
		self.hg18 = open("hg18_clinvar.txt", 'w')
		self.hg18.write("\t".join(new_header) + "\n")
		return self.hg19, self.hg38, self.hg18, self.hg19_bed, self.hg38_bed

	def readin_omim(self):
		with open (self.omim, 'r') as oh:
			for line in oh:
				if line.startswith("#") or re.match(r'\s', line):
					continue
				else:
					terms = line.strip().split("\t")[0].split(",") 
					#If phenotype MIM id is missing in the first field, 
					#it can be found in the third field instead, e.g. 613792
					omim_id_r = re.search(self.omim_re_mm, terms[-1])
					omim_id = omim_id_r.group(1) if omim_id_r else line.strip().split("\t")[2]
					omim_name = ",".join(terms[0:-1]).lower()
					self.omim_dict[omim_name] = omim_id
		return self.omim_dict


	def add_omim_phenotype(self, n):
		if n.lower() in self.omim_dict:
			self.added_names.add(n)
			return self.omim_dict[n.lower()]
		else:
			return self.empty_field

	def parse_clinvar(self):
		self.check_if_corrupted()
		self.initialise_output()
		self.readin_omim()
		file_stream = gzip.open(self.clinvar, 'rt')
		header = file_stream.readline().strip().split("\t")
		for line in file_stream:
			lines = line.strip().split("\t")
			dbsnp_dbvar_id = lines[header.index("RS# (dbSNP)")]
			#Add "rs" prefix if dbSNP identifier present
			if (dbsnp_dbvar_id != "-1"):
				dbsnp_dbvar_id = "rs" + dbsnp_dbvar_id
			#If dbsnp missing, check if dbvar id possibly present
			elif (dbsnp_dbvar_id == "-1" and lines[header.index("nsv/esv (dbVar)")] != "-"):
				dbsnp_dbvar_id = lines[header.index("nsv/esv (dbVar)")]
			else:
				dbsnp_dbvar_id = self.empty_field
			cyto_location = lines[header.index("Cytogenetic")]
			#print(cyto_location)
			if cyto_location == "-":
				cyto_location = self.empty_field
			chromosome = lines[header.index("Chromosome")]
			#do not include if no chromosome given
			if (chromosome == "na" or chromosome == "-1"):
				continue
			#Check if start and end are positive integers.
			start = lines[header.index("Start")]
			end = lines[header.index("Stop")]
			try:
				start_n = int(start)
				end_n = int(end)
				if (start_n > 0 and end_n > 0):
					pass
				else:
					continue
			except:
				continue
			ref_allele = lines[header.index("ReferenceAlleleVCF")]
			alt_allele = lines[header.index("AlternateAlleleVCF")]
			if (ref_allele == "na" or ref_allele == "-"):
				ref_allele = self.empty_field
			if (alt_allele == "na" or alt_allele == "-"):
				alt_allele = self.empty_field
			effect = lines[header.index("ClinicalSignificance")]
			effects = [[x.strip() for x in e.strip().lower().split("/")] for e in effect.split(",")]
			effects_flat = [item for sublist in effects for item in sublist]
			gene = lines[header.index("GeneSymbol")]
			if (self.gene and "subset" in gene):
				continue
			elif (self.gene and "cover" in gene):
				continue
			#Skip non-pathogenic if "Uncertain significance" flag on
			if self.uncertain:
				if ("pathogenic" in effects_flat or "likely pathogenic" in effects_flat or "risk factor" in effects_flat):
					pass
				else:
					continue
			#Can be multiple conditions separated by "|":
			disn = lines[header.index("PhenotypeList")].split("|")
			disi = lines[header.index("PhenotypeIDS")].split("|")
			names_output_final = list()
			omim_output_final = list()
			orpha_output_final = list()
			for (a, b) in zip(disn, disi):
			#can be multiple genes separated by ";"
				disease_name = a.split(";")
				disease_ids = b.split(";")
				names_output = list()
				omim_output = list()
				orpha_output = list()
				for (n, i) in zip(disease_name, disease_ids):
					#Check for presence of omim id
					MyResult_omim = re.search(self.omim_re, i)
					MyResult_orphanet = re.search(self.orphanet_re, i)
					omim_id = MyResult_omim.group(1) if MyResult_omim else self.empty_field
					orphanet_id = MyResult_orphanet.group(1) if MyResult_orphanet else self.empty_field
					#Check if we can assign omim ID ourselves it is missing from the ClinVar file.
					if omim_id == self.empty_field:
						omim_id = self.add_omim_phenotype(n)
					names_output.append(n)
					omim_output.append(omim_id)
					orpha_output.append(orphanet_id)
					names_output1 = ";".join(names_output)	
					omim_output1 = ";".join(omim_output)	
					orpha_output1 = ";".join(orpha_output)
				names_output_final.append(names_output1)
				omim_output_final.append(omim_output1)
				orpha_output_final.append(orpha_output1)
			names_output = "|".join(names_output_final)
			omim_output = "|".join(omim_output_final)
			orpha_output = "|".join(orpha_output_final)
			database_source = "ClinVar"
			VCV = lines[header.index("VariationID")]
			database_link = "http://www.ncbi.nlm.nih.gov/clinvar/variation/" + VCV 
			HGSV_notation = lines[header.index("Name")]
			quality_rating = lines[header.index("ReviewStatus")]
			RCV = lines[header.index("RCVaccession")]
			allele_id = lines[header.index("#AlleleID")]
			assembly = lines[header.index("Assembly")]
			output = [dbsnp_dbvar_id, cyto_location, chromosome, start, end, ref_allele, alt_allele, effect, HGSV_notation, gene, names_output, omim_output, orpha_output, database_source, database_link, VCV, RCV, allele_id, quality_rating] 
			chrom = "chr" + str(chromosome)
			output_bed = [chrom, start, end, ref_allele, alt_allele, dbsnp_dbvar_id, cyto_location, effect, HGSV_notation, gene, names_output, omim_output, orpha_output, database_source, database_link, VCV, RCV, allele_id, quality_rating]
			if assembly == "GRCh37":
				self.hg19.write("\t".join(output) + "\n")
				self.hg19_bed.write("\t".join(output_bed) + "\n")
			elif assembly == "GRCh38":
				self.hg38.write("\t".join(output) + "\n")
				self.hg38_bed.write("\t".join(output_bed) + "\n")
			elif assembly == "NCBI36":
				self.hg18.write("\t".join(output) + "\n")
			else:
				print(assembly)


clinvar_parse = Clinvar(args)
clinvar_parse.parse_clinvar()
