#!/bin/bash
set -e
set -o pipefail
set -u

analysis=$HOME/MendelVar/ontology/inrich_files
scripts=$HOME/bin/mendelvar_standalone/ontology/inrich_files

mkdir -p $analysis
cd $analysis

#DO
cp ../do/do_inrich.txt ./
#DO slim
cp ../do/do_slim_inrich.txt ./

#HPO
cp ../hpo/hpo_inrich.txt ./
#HPO slim
cp ../hpo/hpo_slim_inrich.txt ./

#Freund
cp ../freund/freund_inrich_subset.txt ./

#ConsensusPath
cp ../consensuspath/consensuspath_inrich_subset.txt ./

#GO
cp ../go/go_inrich_subset.txt ./

#GO slim
cp ../go/go_slim_inrich_subset.txt ./

#Pathway Commons
cp ../pathway_commons/pathway_commons_inrich_subset.txt ./

#Reactome
cp ../reactome/reactome_inrich_subset.txt ./