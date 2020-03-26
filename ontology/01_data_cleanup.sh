#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ontology
scripts=$HOME/bin/MendelVar_production/ontology

mkdir -p $analysis && cd $analysis

curl -o "hgnc_complete_set.txt" ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/hgnc_complete_set.txt





