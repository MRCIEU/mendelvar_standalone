#!/usr/bin/env python
from sys import argv
import json
#Create a fake JSON config file, as should be received when running the pipeline on the web.
script, out_f = argv

#Define dictionary with options
options = {
	"assembly" : "GRCh37", #options: GRCh37 | GRCh38 
	"input_type" : "single", #options: single | interval #The single option has to go either with 1k_based_finemapping OR custom_flank
	"interval_generation" : "custom_flank", #options: custom_flank | 1k_based_finemapping (Only for single input_type, otherwise ignored)
	"use_provided_rsids" :"yes", #options yes | no  Option to use rsids provided by user rather than positions in LDlink search. Only for single input and 1k_based_finemapping, otherwise value ignored
	"finemapping_method" : "1k", #options: 1k | 1k_hotspot  #Only for single input and 1k_based_finemapping, otherwise value ignored. 1k_hotspot extends the LD-based interval to the nearest recombination hotspot.
	"LD_metric" : "r2",  #Only for single input and 1k_based_finemapping, otherwise value ignored. #options: r2 | d 
	"LD_threshold" : 0.2,  #Only for single input and 1k_based_finemapping, otherwise value ignored. #options: a decimal from 0 to 1. 
	"flank_left" : 100000, #options: enter a numeric value up to 10 000 000 bp. Only for single input and custom_flank, otherwise ignored
	"flank_right" : 100000, #options: enter a numeric value up to 10 000 000 bp. Only for single input and custom_flank, otherwise ignored
	"target_population" : "EUR", #options: CEU, EUR, AMR, AFR, EAS, SAS 
	"gene_overlap_upstream" : 1000, #options: enter a numeric value up to 20 000 bp. 
	"gene_overlap_downstream" : 1000, #options: enter a numeric value up to 20 000 bp. 
	##Enrichment testing options
	"do" : "yes",  # use Disease ontology or not? Option: yes | no
	"hpo" : "yes",  # use HPO or not? Option: yes | no
	"do_slim" : "yes",  #  use Disease ontology slim or not? Option: yes | no
	"hpo_slim" : "yes",  #  use HPO slim  or not? Option: yes | no
	"freund" : "yes",  # use the Freund set or not? Option: yes | no
	"go" : "no",  # use GO or not? Option: yes | no
	"go_slim" : "no",  # use GO slim  or not? Option: yes | no
	"consensus_path" : "no",  #  use ConsensusPath or not? Option: yes | no
	"pathway_commons" : "no",  #T use PathwayCommons or not? Option: yes | no
	"reactome" : "no",  # use Reactome or not? Option: yes | no
	"inrich_mode" : "gene", #options: interval | gene (ie. target mode) 
	"min_obs_threshold" : 2, #Restrict analysis to the gene sets with at least -z number of genes overlapping with test intervals
	"min_gene_set" : 2, #Restrict analysis to the gene sets with at least (-i) genes 
	"max_gene_set" : 10000  #Restrict analysis to the gene sets with not more than (-j) number of genes

}

out = open(out_f, 'w')
file_dump = json.dump(options, out, indent=4, sort_keys=False)
out.close()

with open (out_f , 'r') as db_h:
	response_dict = json.load(db_h)

print(response_dict)