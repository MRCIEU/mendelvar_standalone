#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd
script, parsed_gff3, apris, output = argv

#Because transcript IDs are not unique
#for loci on X and Y chromosome, cannot
#index dictionary just by transcript ID.
gff3_data = dd(lambda: dd(dict))
transcript_length = dd(int)
with open(parsed_gff3) as gf3:
	for line in gf3:
		lines = line.strip().split()
		position = lines[0] + lines[1] + lines[2]
		#Dictionary indexed by gene_id, transcript_id and position (chrom,start, end)
		gff3_data[lines[6]][lines[5]][position] = lines
		transcript_length[lines[5]] = abs(int(lines[2]) - int(lines[1]))

apris_data = dd(list)
with open(apris) as apris_h:
	for line in apris_h:
		lines = line.strip().split()
		ranking = lines[4]
		if ranking.startswith("PRINCIPAL:"):
			ranking_value = int(ranking.replace("PRINCIPAL:", ""))
			#Check if gene seen before, and if so, if transcript present 
			#in the Gencode data
			if (lines[1] in apris_data and lines[2] in gff3_data[lines[1]]):
				#previous saved hit for the gene
				previous_top = apris_data[lines[1]]
				previous_ranking_value = previous_top[2]
				if previous_ranking_value < ranking_value:
					continue
				elif previous_ranking_value == ranking_value:
					#Compare transcript lengths. If same, pick latest one.
					new_t_length = transcript_length[lines[2]]
					old_t_length = transcript_length[previous_top[1]]
					#If both 0 (i.e. and not found in Gencode database)
					#if (old_t_length == 0 and new_t_length == 0):
					#	continue
					if new_t_length > old_t_length:
						apris_data[lines[1]] = [lines[0], lines[2], ranking_value]
					elif new_t_length == old_t_length:
						apris_data[lines[1]] = [lines[0], lines[2], ranking_value]
					else:
						continue
				else:
					continue
			elif lines[2] in gff3_data[lines[1]]:
				#Save if not seen before
				apris_data[lines[1]] = [lines[0], lines[2], ranking_value]
			else:
				pass
		else:
			continue

#Genes present in ARPPRIS. 
covered_genes = set()
output_fh = open(output, "w")
#Print out the selected transcript ID representing each gene.
for gene in apris_data:
	transcript = apris_data[gene][1]
	if (gene in gff3_data and transcript in gff3_data[gene]):
		gff3 = gff3_data[gene][transcript]
		#Present in both APPRIS and Gencove
		covered_genes.add(gene)
		for entry in gff3:
			to_write = "\t".join(gff3[entry]) #+ "\t" + gene_name
			output_fh.write(to_write + "\n")


#Write out the genes not present in APRIS - select longest transcript, if more than one exists.
additional_transcripts = dict()
for gene in gff3_data:
	if gene in covered_genes:
		continue
	else:
		current_best_length = 0
		for transcript in gff3_data[gene]:
			for position in gff3_data[gene][transcript]:
				my_len = transcript_length[transcript]
				if my_len > current_best_length:
					additional_transcripts[gene] = transcript
					current_best_length = my_len
				else:
					pass
for gene in additional_transcripts:
	transcript = additional_transcripts[gene]
	gff3 = gff3_data[gene][transcript]	
	for entry in gff3:
	#Gene name from APRIS. Compare with entries from Gencode
		to_write = "\t".join(gff3[entry]) #+ "\t" + gene_name
		output_fh.write(to_write + "\n")


output_fh.close()