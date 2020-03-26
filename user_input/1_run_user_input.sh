#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/user_input
scripts=$HOME/bin/mendelvar_standalone/user_input
liftover=$HOME/bin/liftover

mkdir -p $analysis

cd $analysis

##MANUAL CHANGES
#LDlink - requires registering for personal API key. 
#A user can run only one concurrent job using the same LDlink API key
#Change the value of the 'token' variable on line 110 in the script generate_ldlink_interval.py
#to your personal LDlink API key.

#Create a fake JSON file with your run config by modifying the values 
#inside the create_fake_json_config_file.py below
#All the JSON keys have to have a value, 
#but depending on the combination of settings some of them will be ignored. E.g.
#selection of 'interval' input_type will result in ignoring the values of the next 7 fields,
#which become irrelevant.

python $scripts/create_fake_json_config_file.py config.json

#Input file has to be tab-separated, otherwise same format as for MendelVar webserver.
#Run MendelVar:
#1st argument - JSON file with config
#2nd argument - input file

python $scripts/run_mendelvar.py config.json sample_input.bed

#The results will be written to the folder $HOME/MendelVar_user_input/<input_filename> by default.
#Provided inside the directory:
#sample JSON run file: config.json
#sample input file: sample_input.bed