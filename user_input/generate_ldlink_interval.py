#!/usr/bin/env python
from sys import argv
import os, csv, json
import pandas as pd
import numpy as np
import subprocess
import requests, os
from bisect import bisect_left
from collections import defaultdict as dd
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from check_user_input import run_liftover as run_liftover
from check_user_input import liftover_2_final as liftover_2_final
from check_user_input import gtf_2_liftover as gtf_2_liftover
#It does not matter if we specify D or r2' in Ldlink, both are always outputted.

def process_ldlink_output(response, threshold, col_metric, lines, results, logging):
    r = response.content
    my_result = r.decode("ascii").split("\n")
    all_positions = set()
    #Remove empty lines
    my_result = [x for x in my_result if x]
    header = my_result[0].strip().split("\t")
    #If we did not find the correct header indicating that the lookup was successful, terminate and move onto the next position.
    if header != ['RS_Number', 'Coord', 'Alleles', 'MAF', 'Distance', 'Dprime', 'R2', 'Correlated_Alleles', 'FORGEdb', 'RegulomeDB', 'Function']:
        logging.write("WARNING: The following position not found in the 1k reference. Skipping the variant: %s\n" % lines[3])
        return results
    #skipping header
    for entry in my_result[1:]:
        entries = entry.strip().split("\t")
        coordinates = entries[header.index("Coord")].split(":")
        chrom = coordinates[0]
        pos = coordinates[1]
        val_ld = float(entries[header.index(col_metric)])
        if val_ld >= threshold:
            #Collect all the positions returned above a certain LD threshold
            all_positions.add(int(pos))
    #Find the smallest and highest position above the LD threshold
    all_positions_sorted = sorted(all_positions)
    min_position = all_positions_sorted[0]
    max_position = all_positions_sorted[-1]
    results[lines[3]] = (chrom, min_position, max_position)
    return results

def load_hapmap(app_data):
    my_hapmap_file = app_data + "ref_panel/hotspots/all_recombination_hotspots_grch37_sorted.bed"
    hapmap_coordinates = dd(list)
    with open(my_hapmap_file) as mhf:
        for line in mhf:
            lines = line.strip().split()
            hapmap_coordinates[lines[0]].append(int(lines[1]))
    #Convert to a numpy array
    for chrom in hapmap_coordinates:
        my_numpy = np.array(hapmap_coordinates[chrom])
        hapmap_coordinates[chrom] = my_numpy
    return hapmap_coordinates

def get_right_b(max_position, my_numpy_array, max_array):
    if max_position in my_numpy_array:
        right_b = max_position
    elif max_position < max_array:
        right_b = my_numpy_array[my_numpy_array > max_position].min()  
    else:
        right_b = max_position
    return right_b

#Extend to the nearest recombination hotspot (averaged across all pops)
def extend_hotspot(results, my_pars, logging, app_data):
    logging.write("INFO: Extending the LDlink-generated intervals to the nearest non-overlapping recombination hotspots from HapMap3.\n")
    extended_results = {}
    #Read HapMap3 into a dictionary indexed by chromosome and values being a sorted list of positions
    hapmap_coordinates = load_hapmap(app_data)
    for my_id in results:
        (chrom, min_position, max_position) = results[my_id]
        #Find the relevant numpy array:
        if chrom != "chrY":
            my_numpy_array = hapmap_coordinates[chrom]
            min_array = my_numpy_array[0]
            max_array = my_numpy_array[-1]
            max_position = int(max_position)
            min_position = int(min_position)
            #If two positions within a recombination hotspot, do not extend.
            if min_position in my_numpy_array:
                left_b = min_position
                right_b = get_right_b(max_position, my_numpy_array, max_array)
            elif min_position > min_array:
                left_b = my_numpy_array[my_numpy_array < min_position].max()
                right_b = get_right_b(max_position, my_numpy_array, max_array) 
            elif min_position < min_array:
                left_b = min_position
                right_b = get_right_b(max_position, my_numpy_array, max_array)  
            extended_results[my_id] = (chrom, left_b, right_b)
        else:
            extended_results[my_id] = (chrom, min_position, max_position)
    return extended_results

def get_request(url) : 
    i = 0 
    while i < 3 : 
        try: 
            output = requests.get(url, timeout=10000)
            return output  
        except requests.exceptions.RequestException: 
            i += 1 

def requests_ldlink(my_pars, logging, final_input_file):
    #Dictionary indexed by variant ID (the third column) holding the new min and max values for the interval
    results = dict()
    my_pop = my_pars["target_population"]
    token = "my_token"
    my_metric = my_pars["LD_metric"]
    threshold = float(my_pars["LD_threshold"])
    if my_metric == "r2":
        col_metric = "R2"
    else:
        col_metric = "Dprime"
    logging.write("Generating LD-based intervals around subbmitted variants with LDlink\n")
    logging.write("LD metric used: %s, min. LD threshold: %.2f, population: %s\n" % (col_metric, threshold, my_pop))
    
    with open(final_input_file) as fh:
        for line in fh:
            lines = line.strip().split("\t")
            if my_pars["use_provided_rsids"] == "yes":
                query = lines[3].strip()
            else:
                query = lines[0] + ":" + lines[1]
            no_retries = 0
            my_req = "https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=" + query + "&pop=" + my_pop + "&r2_d=" + my_metric + "&token=" + token
            print(my_req)
            response = None
            response = get_request(my_req)
            if response == None:
                logging.write("!!ERROR: We encountered a Timeout error while accessing LDlink. Please try again, or use a custom flank option to generate intervals around each SNP\n")
                raise Exception("We encountered a Timeout error while accessing LDlink. Please try again, or use a custom flank option to generate intervals around each SNP")
            elif (response.status_code != 200):
                r = response.status_code
                print(response.content)
                logging.write("!!ERROR: We encountered a bad status code %s while accessing LDlink. Please try again, or use a custom flank option to generate intervals around each SNP\n"% r)
                raise Exception("We encountered a bad status code %s while accessing LDlink. Please try again, or use a custom flank option to generate intervals around each SNP" % r)
            else:
                results = process_ldlink_output(response, threshold, col_metric, lines, results, logging)
    return results

def create_final_input2(results, logging, final_input_file, my_pars):
    fh = open(final_input_file, 'w')
    for my_id in results:
        fh.write(results[my_id][0] + "\t" + str(results[my_id][1])
        + "\t" + str(results[my_id][2]) + "\t" + my_id + "\n")
    fh.close()
    if my_pars["assembly"] == "GRCh38":
        home = os.environ["HOME"]
        chain_file = str(home) + "/bin/liftover/hg19ToHg38.over.chain"
        logging.write("INFO: Converting the single variants from hg19 to hg38 coordinates with UCSC liftOver tool.\n")
        final_in_file, logging = run_liftover(final_input_file, logging, final_input_file, chain_file)
        num_lines = sum(1 for line in open(final_in_file))
        if num_lines > 0:
            pass
        else:
            logging.write("!!ERROR: Unsuccessful conversion from hg38 to hg19 - no input left. Exiting.\n")
            raise Exception("Unsuccessful conversion from hg38 to hg19 - no input left. Exiting.")
    else:
        fh = open(final_input_file, 'w')
        for my_id in results:
            fh.write(results[my_id][0] + "\t" + str(results[my_id][1])
             + "\t" + str(results[my_id][2]) + "\t" + my_id + "\n")
        fh.close()
    return final_input_file


def run_hotspot_test(my_pars, logging):
    test1 = {"test1": ("chr1", 249077321, 249078149)}
    test1 = extend_hotspot(test1, my_pars, logging, app_data)
    assert test1 == {"test1": ("chr1", 249077321, 249078149)}

    test2 = {"test2": ("chr1", 249077321, 249172191)}
    test2 = extend_hotspot(test2, my_pars, logging, app_data)
    assert test2 == {"test2": ("chr1", 249077321, 249172192)}

    test3 = {"test3": ("chr1", 249077321, 249174684)}
    test3 = extend_hotspot(test3, my_pars, logging, app_data)
    assert test3 == {"test3": ("chr1", 249077321, 249174684)}

    test4 = {"test4": ("chr1", 249077320, 249078149)}
    test4 = extend_hotspot(test4, my_pars, logging, app_data)
    assert test4 == {"test4": ("chr1", 249077287, 249078149)}

    test5 = {"test5": ("chr1", 249077320, 249172191)}
    test5 = extend_hotspot(test5, my_pars, logging, app_data)
    assert test5 == {"test5": ("chr1", 249077287, 249172192)}

    test6 = {"test6": ("chr1", 249077320, 249174684)}
    test6 = extend_hotspot(test6, my_pars, logging, app_data)
    assert test6 == {"test6": ("chr1", 249077287, 249174684)}

    test7 = {"test7": ("chr1", 254990, 249078149)}
    test7 = extend_hotspot(test7, my_pars, logging, app_data)
    assert test7 == {"test7": ("chr1", 254990, 249078149)}

    test8 = {"test8": ("chr1", 254990, 249172191)}
    test8 = extend_hotspot(test8, my_pars, logging, app_data)
    assert test8 == {"test8": ("chr1", 254990, 249172192)}

    test9 = {"test9": ("chr1", 254990, 249174684)}
    test9 = extend_hotspot(test9, my_pars, logging, app_data)
    assert test9 == {"test9": ("chr1", 254990, 249174684)}


def main_ldlink(my_pars, logging, final_input_file, app_data):
    #run_hotspot_test(my_pars, logging)
    print ("Requesting LDlink data")
    results = requests_ldlink(my_pars, logging, final_input_file)
    if my_pars["finemapping_method"] == "1k_hotspot":
        print("Extending to nearest hotspot")
        results = extend_hotspot(results, my_pars, logging, app_data)
    print("Creating final input")
    final_input_file = create_final_input2(results, logging, final_input_file, my_pars)
    return final_input_file
