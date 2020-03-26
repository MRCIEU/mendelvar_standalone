#!/usr/bin/env python
import pandas as pd
from sys import argv

script, popfile, outfile = argv
#Full outer join on the 3 tables produced based on the 3 XML Orphanet files.
filename1 = "all_sites_GRCh38.bed"

f1_df = pd.read_csv(popfile, sep="\t", index_col=False, dtype=str, na_filter = False, header=None, names=["chrom_GChr37", "start_GChr37", "end_GChr37", "rsid"])
f2_df = pd.read_csv(filename1, sep="\t", index_col=False, dtype=str, na_filter = False, header=None, names=["chrom_GChr38", "start_GChr38", "end_GChr38", "rsid"])


df_merge1 = pd.merge(f2_df, f1_df, on='rsid', how='inner')

df_merge1.to_csv(outfile, sep="\t", na_rep="NA", index=False,  encoding='utf-8')

