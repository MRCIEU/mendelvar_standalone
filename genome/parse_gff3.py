#!/usr/bin/env python
from sys import argv
from collections import defaultdict as dd

script, gff3, output_cds, output_transcript = argv

cds_h = open(output_cds, 'w')
t_h = open(output_transcript, 'w')

cds_data = dd(lambda: dd(dict))
with open(gff3) as gff3_h:
	for line in gff3_h:
		if line.startswith("#"):
			pass
		else:
			lines = line.strip().split("\t")
			if (lines[2] == "transcript" or lines[2] == "start_codon" or lines[2] == "stop_codon"):
				my_chrom = lines[0].replace("chr", "")
				my_type = "CDS"
				my_start = lines[3]
				my_end = lines[4]
				my_strand = lines[6]
				my_fields = lines[8].split(";")
				temp_dict = {}
				temp_dict = {x.split("=")[0]: x.split("=")[1] for x in my_fields}
				temp_dict["my_chrom"] = my_chrom
				temp_dict["my_type"] = my_type
				temp_dict["my_strand"] = my_strand
				transcript_id = temp_dict["transcript_id"].split(".")[0]
				gene_id = temp_dict["gene_id"].split(".")[0]
				gene_name = temp_dict["gene_name"].split(".")[0]
				hgnc_id = temp_dict.get("hgnc_id", "NA").replace("HGNC:", "")
				cds_data[transcript_id][my_chrom].update(temp_dict)
				if (lines[2] == "transcript"):
					my_type = "transcript"
					out_line = [my_chrom, my_start, my_end, my_strand, my_type, transcript_id, gene_id, gene_name, hgnc_id]
					t_h.write("\t".join(out_line) + "\n")
				elif (lines[2] == "start_codon" and my_strand == "+"):
					cds_data[transcript_id][my_chrom].update({"my_start": lines[3]})
				elif (lines[2] == "start_codon" and my_strand == "-"):
					cds_data[transcript_id][my_chrom].update({"my_end": lines[4]})

				elif (lines[2] == "stop_codon" and my_strand == "+"):
					cds_data[transcript_id][my_chrom].update({"my_end": lines[4]})
				elif (lines[2] == "stop_codon" and my_strand == "-"):
					cds_data[transcript_id][my_chrom].update({"my_start": lines[3]})
	#print(cds_data)
	for t in cds_data:
		for p in cds_data[t]:
			all_fields = cds_data[t][p] 
			if ("my_start" in all_fields and "my_end" in all_fields):
				out_line = [all_fields["my_chrom"], all_fields["my_start"],
				all_fields["my_end"], all_fields["my_strand"], all_fields["my_type"], 
				all_fields["transcript_id"].split(".")[0], all_fields["gene_id"].split(".")[0],
				all_fields["gene_name"], all_fields.get("hgnc_id", "NA")]
				cds_h.write("\t".join(out_line) + "\n")


cds_h.close()
t_h.close