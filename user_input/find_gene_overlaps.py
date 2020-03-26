#!/usr/bin/env python
import subprocess, os
import pandas as pd
#Create a custom genome file, which extends the gene intervals by values specified by user.

def extend_gene(start, end, strand, upstream_gene_flank, downstream_gene_flank):
    if strand == "-":
        new_start = start - downstream_gene_flank
        new_end = end + upstream_gene_flank
    else:
        new_start = start - upstream_gene_flank
        new_end = end + downstream_gene_flank
    if new_start < 1: new_start = 1
    return new_start, new_end

#Either GRCh38 or GRCh37 depending on input.
def create_custom_genome(my_pars, logging, final_input_file, app_data):
    if my_pars["assembly"] == "GRCh37":
        my_genome = app_data + "genome/gencode.v19.annotation_filtered_transcript_appris_final.txt"
    else:
        my_genome = app_data + "genome/gencode.v31.annotation_filtered_transcript_appris_final.txt"
    upstream_gene_flank = int(my_pars["gene_overlap_upstream"])
    downstream_gene_flank = int(my_pars["gene_overlap_downstream"])
    logging.write("INFO: Extending the GENCODE gene definitions to include %d bp upstream and %d bp downstream of the canononical APPRIS transcript.\n" %(upstream_gene_flank, downstream_gene_flank))
    temp_genome = final_input_file + ".genome"
    temp_genome_fh = open(temp_genome, 'w')
    inrich_genome = final_input_file + ".genome.inrich"
    inrich_genome_fh = open(inrich_genome, 'w')
    #Save the original coordinates of genes with hgnc id - to be used in lookups in disease table
    inrich_dict = dict()
    original_coordinates = dict()
    with open (my_genome) as gh:
        for line in gh:
            lines = line.strip().split()
            chrom = "chr" + lines[0]
            old_start = int(lines[1])
            old_end = int(lines[2])
            strand = lines[3]
            ensg = lines[6]
            name = lines[7]
            hgnc_id = lines[8]
            start, end = extend_gene(old_start, old_end, strand, upstream_gene_flank, downstream_gene_flank)
            out_line = [chrom, str(start), str(end), strand, ensg, name, hgnc_id]
            temp_genome_fh.write("\t".join(out_line) + "\n")
            #Only interested in genes with Mendelian disease in our master table, which all have hgnc_id.
            if hgnc_id != "NA":
                inrich_dict[hgnc_id] = (chrom, start, end)
                original_coordinates[hgnc_id] = (chrom, str(old_start), str(old_end))
    my_keys = sorted(inrich_dict,key=lambda k: (inrich_dict[k][0], inrich_dict[k][1], inrich_dict[k][2]))
    for m in my_keys:
        inrich_genome_fh.write(inrich_dict[m][0].replace("chr", "") + "\t" + str(inrich_dict[m][1]) + "\t" + str(inrich_dict[m][2]) + "\t" + m + "\n")
    temp_genome_fh.close()
    inrich_genome_fh.close()
    return temp_genome, inrich_genome, original_coordinates

def create_giggle_index(logging, final_input_file, temp_genome):
    compressed_file = temp_genome + ".bed" + ".gz"
    cmd = 'cat ' + temp_genome + " | sort --buffer-size 2G -k1,1 -k2,2n -k3,3n | bgzip -c > " + compressed_file
    rcode = subprocess.call(cmd, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem sorting the file with gene coordinates. Exiting.\n")
        raise Exception("We encountered a problem sorting the file with gene coordinates. Exiting.")
    giggle_index = compressed_file + "_index"
    cmd2 = "giggle index -f -i " + compressed_file + " -o " + giggle_index
    rcode = subprocess.call(cmd2, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem indexing the file with gene coordinates. Exiting.\n")
        raise Exception("We encountered a problem indexing the file with gene coordinates. Exiting.")
    return giggle_index

def run_giggle(logging, temp_genome_index, final_input_file, my_pars):
    #Bgzip the query file - required by giggle
    compressed_final = final_input_file + ".gz"
    query_results = final_input_file + "_vs_" + temp_genome_index
    cmd2 = "giggle search -i " + temp_genome_index + " -q " + compressed_final + " -v -o  >" + query_results
    rcode = subprocess.call(cmd2, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem querying the intervals against the gene coordinates. Exiting.\n")
        raise Exception("We encountered a problem querying the intervals against the gene coordinates. Exiting.")
    logging.write("INFO: Successfully queried the user intervals against the gene coordinates from %s \n" % my_pars["assembly"])
    return query_results

def process_giggle_results(query_results, my_pars, inrich_genome, logging):
    upstream_gene_flank = -int(my_pars["gene_overlap_upstream"])
    downstream_gene_flank = -int(my_pars["gene_overlap_downstream"])
    matching_genes = "interval_vs_gene.txt"
    mh = open(matching_genes, 'w')
    header = ["ID", "chrom", "interval_start", "interval_end", "gene_name", "hgnc_id", "Ensembl_gene_id", "gene_start", "gene_end", "strand"]
    mh.write("\t".join(header) + "\n")
    with open(query_results) as giggle_fh:
        for line in giggle_fh:
            if line.startswith("##"):
                lines = line.strip().split("\t")    
                my_rsid = lines[3]
                interval_start = lines[1]
                interval_end = lines[2]
            else:
                lines = line.strip().split("\t")
                strand = lines[3]
                ensg = lines[4]
                name = lines[5]
                hgnc_id = lines[6]
                if hgnc_id != "NA":
                    (chrom, start, end) = inrich_genome[hgnc_id]
                    output_line = [my_rsid, chrom, interval_start, interval_end, name, hgnc_id, ensg, str(start), str(end), strand]
                    mh.write("\t".join(output_line) + "\n")
    mh.close()
    return matching_genes

def check_overlap_empty(query_results, logging):
    if os.stat(query_results).st_size == 0:
        logging.write("!!ERROR:  We found no overlaps between input intervals and genes in the genome. Exiting.\n")
        raise Exception("We found no overlaps between input intervals and genes in the genome. Exiting.")

def check_df_empty(df, logging):
    if df.empty:
        logging.write("!!ERROR: We found no overlaps between input intervals and genes causing Mendelian disease. Exiting.\n")
        raise Exception("We found no overlaps between input intervals and genes causing Mendelian disease. Exiting.")

def merge_with_disease(matching_genes, app_data, logging):
    #Read in a disease table 
    in_disease = app_data + "disease/g_d_integrated_basic_desc_added2_hpo_removed.txt"
    disease_table = pd.read_csv(in_disease, sep="\t", index_col=False, na_filter = False, dtype=object)
    disease_table.drop(columns=["hgnc_gene_name", "ensg_gene_name"], inplace=True, axis=1)
    matched_genes = pd.read_csv(matching_genes, sep="\t", index_col=False, na_filter = False, dtype=object)
    df_merge = pd.merge(matched_genes, disease_table, on='hgnc_id', how='inner')
    check_df_empty(df_merge, logging)
    #Merge with gene synonyms
    hgnc = app_data + "ontology/hgnc_complete_set.txt"
    hgnc_table = pd.read_csv(hgnc, sep="\t", index_col=False, na_filter = False, dtype=object)
    hgnc_table = hgnc_table[['hgnc_id', 'alias_symbol']]
    hgnc_table['hgnc_id'].replace(to_replace=r'^HGNC:', value='', regex=True, inplace=True)
    hgnc_table['alias_symbol'].replace(to_replace=r'^\s*$', value='NA', regex=True, inplace=True)
    filename = "disease_overlap.txt"
    df_merge2 = pd.merge(df_merge, hgnc_table, on="hgnc_id", how="left")
    alt = df_merge2['alias_symbol']
    df_merge2.drop(labels=['alias_symbol'], axis=1,inplace = True)
    df_merge2.insert(5, 'alias_symbol', alt)
    df_merge2['chrom'] = pd.Categorical(df_merge2['chrom'], ["chr1", "chr2", "chr3", "chr4",
    "chr5", "chr6", "chr7", "chr8","chr9", "chr10", "chr11", "chr12", "chr13","chr14", "chr15",
    "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"])
    df_merge2['gene_start'] = df_merge2['gene_start'].astype('int64')
    df_merge2['gene_end'] = df_merge2['gene_end'].astype('int64')
    df_merge2['interval_start'] = df_merge2['interval_start'].astype('int64')
    df_merge2['interval_end'] = df_merge2['interval_end'].astype('int64')
    df_merge2.sort_values(['chrom', 'interval_start', 'interval_end', 'gene_start', "gene_end"], ascending=True, inplace=True)
    df_merge2.to_csv(filename, sep="\t", na_rep="NA", index=False, header=True)
    return filename
    

def run_gene_overlaps(my_pars, logging, final_input_file, app_data):
    print ("Creating custom genome version for INRICH")
    temp_genome, inrich_genome, coordinates_dict = create_custom_genome(my_pars, logging, final_input_file, app_data)
    print ("Creating giggle genome index")
    temp_genome_index = create_giggle_index(logging, final_input_file, temp_genome)
    print ("Running giggle")
    query_results = run_giggle(logging, temp_genome_index, final_input_file, my_pars)
    print ("Checking if gene overlap not empty")
    check_overlap_empty(query_results, logging)
    logging.write("Initialising search for overlaps between input intervals and genes causing Mendelian disease.\n")
    print("Processing giggle results")
    matching_genes = process_giggle_results(query_results, my_pars, coordinates_dict, logging)
    print("Merging with disease data")
    disease_overlap = merge_with_disease(matching_genes, app_data, logging)
    logging.write("Successfully completed search for overlaps between input intervals and genes causing Mendelian disease. " +
    "Found genes written to the disease_overlap.txt file.\n")
    return disease_overlap, inrich_genome, matching_genes