#!/usr/bin/env python
from sys import argv
import os, shutil
from check_user_input import user_init
from check_user_input import check_user_input
from generate_ldlink_interval import main_ldlink
from find_gene_overlaps import run_gene_overlaps
from find_clinvar_overlaps import run_clinvar_overlaps
from run_inrich import run_inrich
from datetime import datetime


startTime = datetime.now()
home = os.environ["HOME"]
user_data = str(home) + "/MendelVar_out/"
app_data = str(home) + "/MendelVar/"
script, my_parameters_file, my_input_file = argv

my_parameters, my_input, logging, user_folder = user_init(my_parameters_file, my_input_file, user_data)

#Read in parameters
my_pars, final_input_file = check_user_input(my_parameters, my_input, app_data, logging)

if my_pars["input_type"] == "single" and my_pars["interval_generation"] == "1k_based_finemapping":
	final_input_file = main_ldlink(my_pars, logging, final_input_file, app_data)

#Return a file with ClinVar variants overlapping our interval. Table #2 - write to file
variant_overlap = run_clinvar_overlaps(final_input_file, my_pars, logging, app_data)
#Table #1 with gene overlaps results - write to file
disease_overlap, inrich_genome, matching_genes = run_gene_overlaps(my_pars, logging, final_input_file, app_data)
#Run INRICH - inrich parsed results and figure
run_inrich(my_pars, logging, final_input_file, inrich_genome, matching_genes, app_data)

print("The execution of the MendelVar pipeline completed in " + str(datetime.now() - startTime))
