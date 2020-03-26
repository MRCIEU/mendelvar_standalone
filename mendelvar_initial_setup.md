## Install GIGGLE, INRICH, LiftOver, bcftools, vcftools in chosen location and add to PATH
[https://github.com/ryanlayer/giggle](https://github.com/ryanlayer/giggle)

[https://atgu.mgh.harvard.edu/inrich/downloads.html](https://atgu.mgh.harvard.edu/inrich/downloads.html)

[https://genome-store.ucsc.edu/](https://genome-store.ucsc.edu)

[http://samtools.github.io/bcftools/bcftools.html](http://samtools.github.io/bcftools/bcftools.html)

[https://vcftools.github.io/](https://vcftools.github.io/)

## Install packages in the conda environment, activate, add to PATH
Needs personalising the `name` and `prefix` in the environment.yml file
```
conda env create environment.yml
```
## Install R dependencies
```
Rscript --vanilla MendelVar_dependencies.R
```

## Prepare hg19 and hg38 gene annotation data
```
sh ./genome/01_genome_cleanup.sh
```

## Prepare 1000 genomes reference panel data (hg19 and hg38) and HapMap3 recombination hotspot data
```
sh ./ref_panel/1_prepare_reference_panel.sh
```

## Populate OMIM data
User needs to be a registered OMIM user with a personal API Key.

apiKey needs to be changed to personal API key in the script `01_omim_build.sh` (lines 12,13,14) and in the script `omim_api.py` (line 21).
```
sh ./omim/01_omim_build.sh
```

## Get HGNC reference for ontology analysis
```
sh ./ontology/01_data_cleanup.sh
```

## Get Freund et al. gene sets
```
sh ./ontology/freund/1_freund_ontology_cleanup.sh
```

## Run the regular update script to get other datasets
First, read through the script and make the required manual changes to it.
```
sh ./user_input/2_mendelvar_regular_update.sh
```