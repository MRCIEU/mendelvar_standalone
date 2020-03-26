#!/usr/bin/env python
import subprocess, re, os
from inrich_wrapper import Inrich
from inrich_wrapper import GeneMap
from inrich_wrapper import GeneSet
from inrich_wrapper import MatchingGenes

def create_inrich_interval(final_input_file):
    inrich_interval = final_input_file + ".inrich"
    ii_fh = open(inrich_interval, "w")
    with open(final_input_file) as fh:
        for line in fh:
            lines = line.strip().split()
            lines[0] = lines[0].replace("chr", "")
            ii_fh.write("\t".join(lines[0:3]) + "\n")
    ii_fh.close()
    return inrich_interval

def write_parsed_inrich(inrich_in, genemap_parsed, gene_set, my_interval_genes_parsed, out, min_obs_threshold):
    with open(out, "w") as output:
        output.write("Gene_set_size" + "\t" + "Overlap_size" + "\t" + "Overlapping_genes_hgnc_id"
        + "\t" +  "Overlapping_genes_hgnc_names" + "\t" + "Empirical_p_value"+ "\t" + "Bootstrapped_p_value" 
         + "\t" + "Gene_set_ID" + "\t" + "Gene_set_description" + "\n")
        #Sort INRICH output by empirical p-value
        my_sorted = sorted(inrich_in.my_data.items(), key=lambda x: (float(x[1][2]), int(x[1][1])))
        for key, values in my_sorted:
            #Skip zero gene overlaps with term
            if int(values[1]) >= int(min_obs_threshold):
                gene_ids = list(gene_set[key][1])
                gene_names = [my_interval_genes_parsed[g] for g in gene_ids if g in my_interval_genes_parsed]
                gene_ids2 = [g for g in gene_ids if g in my_interval_genes_parsed]
                #fdr = fdr_pvals.get(float(values[2]) ,"NA")
                my_line = values[0] + "\t" + values[1] + "\t" + ", ".join(gene_ids2) + "\t" + ", ".join(gene_names) + "\t" + values[2] + "\t" + values[3] + "\t" + key + "\t" + gene_set[key][0] + "\n"
                #if int(values[1]) != len(gene_names):
                #    print (my_line)
                output.write(my_line)
    return out

def parse_inrich_results(inrich_output, ont_source, app_data, matching_genes, my_pars):
    inrich_output = inrich_output + ".out.inrich" 
    hgnc = app_data + "ontology/hgnc_complete_set.txt"
    my_inrich = Inrich(inrich_output)
    inrich_parsed = my_inrich.read_in_inrich()
    #fdr_pvals = inrich_parsed.pvalue()
    my_genemap = GeneMap(hgnc)
    genemap_parsed = my_genemap.read_in_gene_map()
    my_interval_genes = MatchingGenes(matching_genes)
    my_interval_genes_parsed = my_interval_genes.read_in_genes()
    my_geneset = GeneSet(ont_source)
    geneset_parsed = my_geneset.read_in_gene_set()
    out = inrich_output + ".parsed"
    min_obs_threshold = my_pars.get("min_obs_threshold", 2)
    inrich_results = write_parsed_inrich(inrich_parsed, genemap_parsed, geneset_parsed, my_interval_genes_parsed, out, min_obs_threshold)
    return inrich_results

def plot_inrich_results(out, app_data, logging):
    plotr = app_data + "ontology/mendelvar_figure.R"
    print(out)
    cmd = "Rscript --vanilla " + plotr + " " + out
    inrich_figure_png = out.replace(".out.inrich.parsed", ".png")
    inrich_figure_pdf = out.replace(".out.inrich.parsed", ".pdf")
    rcode = subprocess.call(cmd, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem plotting the enrichment results. Exiting.\n")
    else:
        logging.write("INFO: Successfully completed plotting results of enrichment testing. Results written to file: " + inrich_figure_png + "\n")
        logging.write("INFO: Successfully completed plotting results of enrichment testing. Results written to file: " + inrich_figure_pdf + "\n")
    return inrich_figure_png, inrich_figure_pdf

def single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars):
    print ("Running INRICH")
    my_ont_file = os.path.basename(ont_source)
    inrich_output = re.sub('_inrich.*', "", my_ont_file)
    min_obs_threshold = my_pars.get("min_obs_threshold", "2")
    min_gene_set = my_pars.get("min_gene_set", "2")
    max_gene_set = my_pars.get("max_gene_set", "10000")

    logging.write("INFO: Running INRICH using the following ontology: %s\n" % inrich_output)
    cmd = "inrich -c -q 10000 -r 50000 -p 1 -i " + str(min_gene_set) + " -j " +  str(max_gene_set) +  " -z " + str(min_obs_threshold) + " -a " + inrich_interval + " -m " + my_snp_map + " -g " + inrich_genome + " -t " + ont_source + " -o " + inrich_output + mode + "\n"
    rcode = subprocess.call(cmd, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem running ontology enrichment with INRICH. Exiting.\n")
        raise Exception("We encountered a problem running ontology enrichment with INRICH. Exiting.")
    print ("Parsing INRICH results")
    inrich_results = parse_inrich_results(inrich_output, ont_source, app_data, matching_genes, my_pars)
    logging.write("INFO: Successfully completed enrichment testing. Results written to file: " + inrich_results + "\n")
    print ("Plotting INRICH results")
    inrich_figure_png, inrich_figure_pdf = plot_inrich_results(inrich_results, app_data, logging)
    return inrich_results, inrich_figure_png, inrich_figure_pdf

def run_inrich(my_pars, logging, final_input_file, inrich_genome, matching_genes, app_data):
    #Create INRICH interval file
    print ("Creating INRICH intervals")
    inrich_interval = create_inrich_interval(final_input_file)
    #Pick the right SNP map file.
    my_pop = my_pars["target_population"]
    if my_pars["assembly"] == "GRCh38":
        my_snp_map = app_data + "ref_panel/" + my_pop + "_sites_1k_GRCh38.inrich"
    else:
        my_snp_map = app_data + "ref_panel/" + my_pop + "_sites_1k_GRCh37.inrich"
    if my_pars["inrich_mode"] == "gene":
        mode = " -2"
    else:
        mode = " "
    #Go over all possible ontologies

    if my_pars.get("do", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "do_inrich.txt"
        inrich_do, do_figure_png, do_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("do_slim", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "do_slim_inrich.txt"
        inrich_do_slim, do_slim_figure_png, do_slim_figure_pdf  = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("hpo", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "hpo_inrich.txt"
        inrich_hpo, hpo_figure_png, hpo_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
  
    if my_pars.get("hpo_slim", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "hpo_slim_inrich.txt"
        inrich_hpo_slim, hpo_slim_figure_png, hpo_slim_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
  
    if my_pars.get("freund", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "freund_inrich_subset.txt"
        inrich_freund, freund_figure_png, freund_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("go", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "go_inrich_subset.txt"
        inrich_go, go_figure_png, go_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("go_slim", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "go_slim_inrich_subset.txt"
        inrich_go_slim, go_slim_figure_png, go_slim_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("consensus_path", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "consensuspath_inrich_subset.txt"
        inrich_consensus_path, consensuspath_figure_png, consensuspath_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("pathway_commons", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "pathway_commons_inrich_subset.txt"
        inrich_pathway_commons, pathway_commons_figure_png, pathway_commons_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)
    if my_pars.get("reactome", "no") == "yes":
        ont_source = app_data + "ontology/inrich_files/" + "reactome_inrich_subset.txt"
        inrich_reactome, reactome_figure_png, reactome_figure_pdf = single_inrich(inrich_interval, my_snp_map, inrich_genome, ont_source, mode, app_data, matching_genes, logging, my_pars)

