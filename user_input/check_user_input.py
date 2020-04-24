#!/usr/bin/env python
from sys import argv
import os, csv, json
import pandas as pd
import numpy as np
import subprocess
from pathlib import Path

def user_init(my_parameters_file, my_input_file, user_data):
    pd.set_option("mode.chained_assignment", None)
    user_folder =  user_data + os.path.basename(my_input_file) + "/"
    cmd_1 = "mkdir -p " + user_folder
    rcode = subprocess.call(cmd_1, shell=True)
    if rcode == 1:
        raise Exception("We encountered a problem creating a working directory. Exiting.")
    cmd_2 = 'cp ' + my_input_file + " " + user_folder + " && cp " + my_parameters_file + " " + user_folder
    rcode = subprocess.call(cmd_2, shell=True)
    if rcode == 1:
        raise Exception("We encountered a problem while copying the input file. Exiting.")
    os.chdir(user_folder)
    logging = open(user_folder + "mendelvar.log", 'w')
    return os.path.basename(my_parameters_file), os.path.basename(my_input_file), logging, user_folder


def user_init_celery(my_parameters_file, my_input_file, job_id, user_data):
    pd.set_option("mode.chained_assignment", None)
    user_folder =  user_data + job_id + "/"
    os.chdir(user_folder)
    logging = open(user_folder + "mendelvar.log", 'w')
    return os.path.basename(my_parameters_file), os.path.basename(my_input_file), logging, user_folder

def read_in_params(my_parameters, logging):
    with open (my_parameters, 'r') as db_h:
        my_pars = json.load(db_h)
    logging.write("********MendelVar logfile********" + "\n")
    logging.write("!!!!Running MendelVar with the following settings:--->" + "\n")
    my_pars_short = my_pars.copy()
    if my_pars["input_type"] == "interval":
        remove_keys = ["interval_generation", "use_provided_rsids", "finemapping_method", "LD_metric", "LD_threshold", "flank_left", "flank_right"]
        for k in remove_keys:
            my_pars_short.pop(k, None)
    elif my_pars['input_type'] == "single":
        if my_pars['interval_generation'] == "custom_flank":
            remove_keys = ["use_provided_rsids", "finemapping_method", "LD_metric", "LD_threshold"]
            for k in remove_keys:
                my_pars_short.pop(k, None)
        elif my_pars['interval_generation'] == "1k_based_finemapping":
            remove_keys = ["flank_left", "flank_right"]
            for k in remove_keys:
                my_pars_short.pop(k, None)
    my_pars_short.pop("json", None)
    for k, v in my_pars_short.items():
        logging.write("\t" + str(k) + " : " + str(v) + "\n")
    logging.write("!!!!End of settings file--->\n\n")
    return my_pars

def check_number_cols(my_input, logging):
    rows_no = []
    try:
        reader = csv.reader(open(my_input, 'r'),  delimiter='\t', skipinitialspace=True)
        for row in reader:
            #Check for empty lines
            if any(row):
                rows_no.append(len(row))
    except UnicodeDecodeError:
        logging.write("!!ERROR: Invalid, non-text binary input file. Exiting.\n")
        raise Exception("Invalid, non-text binary input file. Exiting.")     
    rows_no_uniq = set(rows_no)
    col_n = list(rows_no_uniq)
    if len(rows_no_uniq) == 0:
        logging.write("!!ERROR: Detected empty input file. Exiting.\n")
        raise Exception("Detected empty input file. Exiting.")
    elif len(rows_no_uniq) > 1:
        rows_no_uniq_concat = ", ".join([str(r) for r in rows_no_uniq])
        logging.write("!!ERROR: Lines with variable number of columns detected in the input: %s columns. Exiting.\n" % rows_no_uniq_concat)
        raise Exception("Lines with variable number of columns detected in the input: %s columns. Exiting." % rows_no_uniq_concat)
    #Only allowing 2 columns to 4 columns
    elif (col_n[0] < 2 or col_n[0] > 4):
        logging.write("!!ERROR: Incorrect number of columns in input: %d - only between 2 and 4 allowed. Exiting.\n" % col_n[0])
        raise Exception("Incorrect number of columns in input: %d - only between 2 and 4 allowed. Exiting." % col_n[0])
    else:
        return col_n[0]

def read_in_df(my_input, names):
    df = pd.read_csv(my_input, sep="\t", index_col=False, na_filter = False, header=None, names=names, dtype=object)
    return df

def read_in_file(my_input, col_n, my_pars, logging):
    if col_n == 2:
        if my_pars["input_type"] == "single":
            names = ["chrom", "pos"]
            df = read_in_df(my_input, names)
        else:
            logging.write("!!ERROR: Incorrect input_type option given: interval - did you mean lead SNPs? Exiting.\n")
            raise Exception("Incorrect input_type option given: interval - did you mean lead SNPs? Exiting.")   
    if col_n == 3:
        if my_pars["input_type"] == "single":   
            names = ["chrom", "pos", "custom_id"]
            df = read_in_df(my_input, names)
        if my_pars["input_type"] == "interval":   
            names = ["chrom", "start", "end"]
            df = read_in_df(my_input, names)      
    if col_n == 4:
        if my_pars["input_type"] == "interval":
            names = ["chrom", "start", "end", "custom_id"]
            df = read_in_df(my_input, names)
        else:
            logging.write("!!ERROR: Incorrect input_type option given: lead SNPs - did you mean interval? Exiting.\n")
            raise Exception("Incorrect input_type option given: lead SNPs - did you mean interval? Exiting.")     
    #Drop empty lines
    df.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
    #Drop duplicate lines
    df.drop_duplicates(inplace=True)
    return df

def check_max(df, logging):
    number_records = len(df.index) 
    if number_records > 10000:
        logging.write("!!ERROR: Number of entries in the input is %d, higher than the max limit of 10000 allowed. Exiting.\n" % number_records)
        raise Exception("Number of entries in the input is %d, higher than the max limit of 10000 allowed. Exiting." % number_records)

def check_max_interval(df, logging):
    if ('start' in df and 'end' in df):
        df["difference"] = abs(df["end"] - df["start"])
        if (df["difference"] > 20000000).any():
            logging.write("!!ERROR: Input intervals exceed the maximum length of 20 Mbp. Exiting.\n")
            raise Exception("Input intervals exceed the maximum length of 20 Mbp. Exiting.")
        df.drop(columns=["difference"], inplace=True)
    return df

def generate_flanks(df, my_pars, logging):
    if (my_pars["input_type"] == "single" and my_pars["interval_generation"] == "custom_flank"):
        df['pos'] = df['pos'].astype('int64')
        start_series = df['pos'] - int(my_pars['flank_left'])
        end_series = df['pos'] + int(my_pars['flank_right'])
        logging.write("INFO: Adding custom flanks around input variants.\n")
        logging.write("INFO: Adding %d bp upstream of the variant\n" % int(my_pars['flank_left']))
        logging.write("INFO: Adding %d bp downstream of the variant\n" % int(my_pars['flank_right']))
        df.drop(columns=["pos"], inplace=True, axis=1)
        df.insert(1, 'start', start_series, allow_duplicates = False)
        df.insert(2, 'end', end_series, allow_duplicates = False)
        #Coerce to 1 if negative value after generating upstream flank
        df.loc[(df.start < 1),'start'] = 1
        if len(df.index) == 0:
            logging.write("!!ERROR: No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.\n")
            raise Exception("No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.")  
        else:
            return df      
    else:
        return df

#Take care of columns data types
def convert_columns(df, logging):
    #First column
    allowed_chroms = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
    if 'chrom' in df.columns:
        df['chrom'].replace(to_replace=r'^chr', value='', regex=True, inplace=True)
        #Check that only 1-22 and X, Y chrom
        to_keep = df['chrom'].isin(allowed_chroms)
        df = df.loc[to_keep]    
        df['chrom'].replace(to_replace=r'^', value='chr', regex=True, inplace=True)
    #Second column
    if "pos" in df.columns:
        #signed integer, max 4294967295
        df[["pos"]] = df[["pos"]].apply(pd.to_numeric, downcast="integer", errors="coerce")
        #Manually drop postions smaller than 1 as when using unsigned coercion they are coerced to big integers.     
        df.loc[(df.pos < 1),'pos'] = np.nan
        df.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        df['pos'] = df['pos'].astype('uint64')
    if "start" in df.columns:
        #signed integer, max 4294967295
        df[["start"]] = df[["start"]].apply(pd.to_numeric, downcast="integer", errors="coerce")
        #Manually drop postions smaller than 1 as when using unsigned coercion they are coerced to big integers.     
        df.loc[(df.start < 1),'start'] = np.nan
        df.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        df['start'] = df['start'].astype('uint64')
    if "end" in df.columns:
        #signed integer, max 4294967295
        df[["end"]] = df[["end"]].apply(pd.to_numeric, downcast="integer", errors="coerce")
        df['end'] = np.where(df['end'] < 1, np.nan, df['end'])
        df.dropna(axis=0, how='any', thresh=None, subset=None, inplace=True)
        df['end'] = df['end'].astype('uint64')
    if len(df.index) == 0:
        logging.write("!!ERROR: No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.\n")
        raise Exception("No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.")          
    return df

#Check if start coordinate smaller than end coordinate
def check_values(df, logging):
    if ('start' in df and 'end' in df):
        df["comparison"] = np.where(df['end'] >= df['start'], df['start'], np.nan)
        df.dropna(inplace=True)
        df.drop(columns=["comparison"], inplace=True)
    if "custom_id" in df:
        unique_ids = df["custom_id"].nunique()
        if len(df.index) != unique_ids:
            logging.write("!!ERROR: Non unique IDs submitted by the user found. Exiting.\n")
            raise Exception("Non unique IDs submitted by the user found. Exiting.")  
    if len(df.index) == 0:
        logging.write("!!ERROR: No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.\n")
        raise Exception("No entries left after filtering on column values. Check if they meet the description in the tutorial. Exiting.")           
    return df

def append_ids(df):
    if "custom_id" in df:
        pass
    else:
        if "start" in df:
            df["custom_id"] = df["chrom"].astype(str) + "_" + df["start"].astype(str) + "_" + df["end"].astype(str)
        else:
            df["custom_id"] = df["chrom"].astype(str) + "_" + df["pos"].astype(str)
    return df

def write_to_file(df, my_input):
    filename = os.path.basename(my_input) + ".sorted"
    df.to_csv(filename, sep="\t", na_rep="NA", index=False, header=False)
    return filename

def sort_columns(df, my_input, logging):
    if ('chrom' in df and 'start' in df and 'end' in df): 
        df.sort_values(['chrom', 'start', 'end'], ascending=True, inplace=True)
        return write_to_file(df, my_input)
    elif ('chrom' in df and 'pos' in df):
        df.sort_values(['chrom', 'pos'], ascending=True, inplace=True)
        return write_to_file(df, my_input)
    else:
        logging.write("!!ERROR: Something went wrong. Malformed input for sorting. Exiting.\n")
        raise Exception("Something went wrong. Malformed input for sorting. Exiting.")

def check_empty(query_results, logging):
    if os.stat(query_results).st_size == 0:
        logging.write("!!ERROR: We found no matching positions when transferring input variants from GRCh38 to GrCh37 coordinates. Exiting.\n")
        raise Exception("We found no matching positions when transferring input variants from GRCh38 to GrCh37 coordinates. Exiting.")
    else:
        return 1

def gtf_2_liftover(sorted_file):
    outfile = sorted_file + ".in.bed"
    outfile_h = open(outfile, 'w')
    with open(sorted_file, 'r') as sof:
        for line in sof:
            if line:
                lines = line.strip().split("\t")
                if len(lines) == 3:
                    lines.append(lines[2])
                    lines[2] = str(int(lines[1]))
                    lines[1] = str(int(lines[1]) - 1)
                else:
                    lines[1] = str(int(lines[1]) - 1)
                outfile_h.write(" ".join(lines) + "\n")
        outfile_h.close()
    return outfile

def gtf_2_giggle(sorted_file):
    outfile = sorted_file + ".in.bed"
    outfile_h = open(outfile, 'w')
    with open(sorted_file, 'r') as sof:
        for line in sof:
            if line:
                lines = line.strip().split("\t")
                if len(lines) == 3:
                    lines.append(lines[2])
                    lines[2] = str(int(lines[1]))
                    lines[1] = str(int(lines[1]))
                outfile_h.write("\t".join(lines) + "\n")
        outfile_h.close()
    return outfile
    
def liftover_2_final(mapped, final_in_file):
    outfile_h = open(final_in_file, 'w')
    with open(mapped) as ifh:
        for line in ifh:
            if line:
                lines = line.strip().split("\t")
                lines[1] = str(int(lines[1]) + 1)
                outfile_h.write("\t".join(lines) + "\n")
    outfile_h.close()


def run_liftover(sorted_file, logging, final_in_file, chain_file):
        mapped = final_in_file + ".mapped"
        unmapped = final_in_file + ".unmapped"
        infile = gtf_2_liftover(sorted_file)
        cmd = 'liftOver ' + infile + " " + chain_file + " " + mapped + " " + unmapped
        rcode = subprocess.call(cmd, shell=True)
        if rcode == 0:
            check_empty(mapped, logging)
            liftover_2_final(mapped, final_in_file)
            unmapped_h = open(unmapped)
            unmapped_all = unmapped_h.read()
            logging.write("*** Start of unmapped positions.\n")  
            logging.write(unmapped_all)
            unmapped_h.close()
            logging.write("*** End of unmapped positions.\n")  
            return final_in_file, logging     
        else:
            logging.write("We encountered a problem converting from hg38 to hg19 coordinates. Exiting.\n")
            raise Exception("We encountered a problem converting from hg38 to hg19 coordinates. Exiting.")

def parse_1k_giggle(query_results, final_in_file, logging):
    vh = open(final_in_file, 'w')
    with open(query_results) as giggle_fh:
        for line in giggle_fh:
            if line.startswith("##"):
                lines = line.strip().split("\t")    
                ID = lines[3]
                old_start = lines[1]
                old_end = lines[2]
            else:
                lines = line.strip().split("\t")
                chrom = lines[0]
                new_start = lines[5]
                new_end = lines[6]
                vh.write(chrom + "\t" + new_start + "\t" + new_end + "\t" + ID + "\n")
    vh.close()
    return final_in_file

def lookup_1k_GRCh38(sorted_file, logging, my_pars, app_data):
    #Find the right population index
    my_pop = my_pars["target_population"]
    pop_index = app_data + "ref_panel/" + my_pop + "_sites_GRCh38.bed_index"
    my_bed = pop_index.replace("_index", ".gz")
    #Create a symbolic link to the 1k bed file required by giggle
    if not os.path.exists(os.path.basename(my_bed)):
        create_link = "ln -s " + my_bed + " " + os.path.basename(my_bed)
        rcode = subprocess.call(create_link, shell=True)
        if rcode == 1:
            logging.write("!!ERROR: We encountered a problem linking the 1k variant index. Exiting.\n")
            raise Exception("We encountered a problem linking the 1k variant index. Exiting.")
    query_file = gtf_2_giggle(sorted_file)
    compressed_final = query_file + ".gz"
    cmd = "cat " + query_file + " | bgzip -c > " + compressed_final
    rcode = subprocess.call(cmd, shell=True)
    if rcode == 1:
        logging.write("!!ERROR: We encountered a problem compressing the query file. Exiting.\n")
        raise Exception("We encountered a problem compressing the query file. Exiting.")
    query_results = sorted_file + "_vs_1k_" + my_pop
    cmd2 = "giggle search -i " + pop_index + " -q " + compressed_final + " -v -o  >" + query_results
    rcode = subprocess.call(cmd2, shell=True)
    if rcode == 1:
        logging.write("!!Error: We encountered a problem querying the input variants against the 1k variant coordinates. Exiting.\n")
        raise Exception("We encountered a problem querying the input variants against the 1k variant coordinates. Exiting.")
    else:
        check_empty(query_results, logging)
    return query_results

def create_final_input(sorted_file, my_pars, app_data, logging):
    final_in_file = sorted_file + ".bed"
    if (my_pars['assembly'] == "GRCh38" and my_pars['input_type'] == "single" and my_pars['interval_generation'] == "1k_based_finemapping" and my_pars["use_provided_rsids"] == "no"):
        home = os.environ["HOME"]
        chain_file = str(home) + "/bin/liftover/hg38ToHg19.over.chain"
        logging.write("INFO: Converting the single variants from hg38 to hg19 coordinates with UCSC liftOver tool.\n")
        final_in_file, logging = run_liftover(sorted_file, logging, final_in_file, chain_file)
        num_lines = sum(1 for line in open(final_in_file))
        if num_lines > 1:
            pass
        else:
            logging.write("!!ERROR: Unsuccessful conversion from hg39 to hg19 - no input left. Exiting.\n")
            raise Exception("Unsuccessful conversion from hg39 to hg19 - no input left. Exiting.")    
    elif (my_pars['input_type'] == "single" and my_pars['interval_generation'] == "1k_based_finemapping"):
        fh_out = open(final_in_file, 'w')
        with open(sorted_file) as sh:
            for line in sh:
                lines = line.strip().split()
                lines.append(lines[2])
                lines[2] = lines[1]
                fh_out.write("\t".join(lines) + "\n")
        fh_out.close()
    else:
        cmd = 'cp ' + sorted_file + " " + final_in_file
        rcode = subprocess.call(cmd, shell=True)
        if rcode == 1:
            logging.write("!!ERROR: We encountered a problem creating a final input file. Exiting.\n")
            raise Exception("We encountered a problem creating a final input file. Exiting.")
    return final_in_file


def check_user_input(my_parameters, my_input, app_data, logging):
    my_pars = read_in_params(my_parameters, logging)
    print("Checking number of columns")
    col_n = check_number_cols(my_input, logging)
    print("Loading file")
    df = read_in_file(my_input, col_n, my_pars, logging)
    print("Checking number of entries")
    check_max(df, logging)
    print("Converting columns")
    df = convert_columns(df, logging)
    print("Generating intervals")
    df = generate_flanks(df, my_pars, logging)
    print("Checking interval length")
    check_max_interval(df, logging)
    print("Checking correctness of interval values")
    df = check_values(df, logging)
    print("Appending IDs")
    df = append_ids(df)
    print("Sorting columns")
    sorted_file = sort_columns(df, my_input, logging)
    #Convert to GrCh37 when needed
    #Make sure 4 columns
    print("Creating final input")
    final_input = create_final_input(sorted_file, my_pars, app_data, logging)
    return my_pars, final_input

