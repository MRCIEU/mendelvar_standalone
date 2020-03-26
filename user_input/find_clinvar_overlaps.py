#!/usr/bin/env python
import subprocess, os

def run_giggle(final_input_file, my_pars, logging, app_data):
    print ("Running giggle")
    logging.write("INFO: Initialising search for overlaps between input intervals and ClinVar variants.\n")
    if my_pars["assembly"] == "GRCh37":
        my_index = "hg19_clinvar.bed_index"
        my_variant_index = app_data + "hg19_clinvar.bed_index"
    else:
        my_index = "hg38_clinvar.bed_index"
    
    my_variant_index = app_data + "clinvar/" + my_index
    my_bed = app_data + "clinvar/" + my_index.replace("_index", ".gz")
    #Create a symbolic link to the ClinVar bed file required by giggle
    if not os.path.exists(os.path.basename(my_bed)):
        create_link = "ln -s " + my_bed + " " + os.path.basename(my_bed)
        rcode = subprocess.call(create_link, shell=True)
        if rcode == 1:
            logging.write("!!ERROR: We encountered a problem linking the ClinVar variant index. Exiting.\n")
            raise Exception("We encountered a problem linking the ClinVar variant index. Exiting.")
    compressed_final = final_input_file + ".gz"
    cmd = "cat " + final_input_file + " | bgzip -c > " + compressed_final
    rcode = subprocess.call(cmd, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem compressing the query file. Exiting.\n")
        raise Exception("We encountered a problem compressing the query file. Exiting.")
    query_results = final_input_file + "_vs_" + my_index
    cmd2 = "giggle search -i " + my_variant_index + " -q " + compressed_final + " -v -o  >" + query_results
    rcode = subprocess.call(cmd2, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem querying the intervals against the ClinVar variant coordinates. Exiting.\n")
        raise Exception("We encountered a problem querying the intervals against the ClinVar variant coordinates. Exiting.")
    logging.write("INFO: Successfully queried the user intervals against the ClinVar variant coordinates from %s. Found variants written to the variant_overlap.txt file\n" % my_pars["assembly"])
    return query_results

def check_overlap_empty(query_results, logging):
    print("Checking if overlaps not empty")
    if os.stat(query_results).st_size == 0:
        logging.write("INFO: We found no overlaps between input intervals and variants in ClinVar.\n")
        return 0
    else:
        return 1

def initiate_output_file():
    print ("Initiating output file")
    out_file = "variant_overlap.txt"
    vh = open(out_file, 'w')
    header = ["ID", "chrom", "interval_start", "interval_end", "dbsnp_dbvar_id", "variant_start", "variant_end", "cyto_location", 
        "ref_allele", "alt_allele", "effect", "HGSV_notation", "gene", "disease_name(s)", "disease_omim", "disease_orphanet",
        "clinvar_link", "VCV", "RCV", "allele_id", "quality_rating"]
    vh.write("\t". join(header) + "\n")
    vh.close()
    return out_file

def process_giggle_results(query_results, logging, variant_overlap):
    print ("Processing Giggle results")
    all_lines = list()
    chrom_order = {"chr1":1, "chr2":2, "chr3":3, "chr4":4,
    "chr5":5, "chr6":6, "chr7":7, "chr8":8,"chr9":9, "chr10":10, "chr11":11, "chr12":12, "chr13":13,"chr14":14, "chr15":15,
    "chr16":16, "chr17":17, "chr18":18, "chr19":19, "chr20":20, "chr21":21, "chr22":22, "chrX":23, "chrY":24}
    vh = open(variant_overlap, 'a')
    with open(query_results) as giggle_fh:
        for line in giggle_fh:
            if line.startswith("##"):
                lines = line.strip().split("\t")    
                ID = lines[3]
                interval_start = lines[1]
                interval_end = lines[2]
            else:
                lines = line.strip().split("\t")
                chrom = lines[0]
                dbsnp_dbvar_id = lines[5]
                variant_start = lines[1]
                variant_end = lines[2]
                cyto_location = lines[6]
                ref_allele = lines[3]
                alt_allele = lines[4]
                effect = lines[7]
                HGSV_notation = lines[8]
                gene = lines[9]
                disease_name = lines[10]
                disease_omim = lines[11]
                disease_orphanet = lines[12].replace("ORPHA", "")
                clinvar_link = lines[14]
                VCV = lines[15]
                RCV = lines[16].replace("RCV","")
                allele_id = lines[17]
                quality_rating = lines[18]
                outline = (ID, chrom, int(interval_start), int(interval_end), dbsnp_dbvar_id, int(variant_start), int(variant_end), cyto_location, 
                ref_allele, alt_allele, effect, HGSV_notation, gene, disease_name, disease_omim, disease_orphanet,
                clinvar_link, VCV, RCV, allele_id, quality_rating)
                all_lines.append(outline)
    all_lines2 = sorted(all_lines, key = lambda x: (chrom_order[x[1]], x[2], x[3], x[5], x[6]))
    for a in all_lines2:
        a = [str(x) for x in a]
        vh.write("\t".join(a) + "\n")
    vh.close()
    logging.write("INFO: Successfully completed search for overlaps between input intervals and ClinVar variants. " +
    "Found variants written to the variant_overlap.txt file\n")
    return variant_overlap
    

def run_clinvar_overlaps(final_input_file, my_pars, logging, app_data):
    query_results = run_giggle(final_input_file, my_pars, logging, app_data)
    present_overlap = check_overlap_empty(query_results, logging)
    variant_overlap = initiate_output_file()
    if present_overlap:
        variant_overlap = process_giggle_results(query_results, logging, variant_overlap)
    return variant_overlap

