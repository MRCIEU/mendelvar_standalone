#!/usr/bin/env python
import pandas as pd
#Full outer join on the 3 tables produced based on the 3 XML Orphanet files.
filename1 = "orphanet_xml1_parsed"
filename4 = "orphanet_xml4_parsed"
filename6 = "orphanet_xml6_parsed"
f1_df = pd.read_csv(filename1, sep="\t", index_col=False, dtype=str, na_filter = False, encoding="latin")
f4_df = pd.read_csv(filename4, sep="\t", index_col=False, dtype=str, na_filter = False, encoding="latin")
f6_df = pd.read_csv(filename6, sep="\t", index_col=False, dtype=str, na_filter = False, encoding="latin")

df_merge1 = pd.merge(f1_df, f6_df, on='orphanet_id', how='outer')
df_merge2 = pd.merge(df_merge1, f4_df, on='orphanet_id', how='outer')
df_merge2.to_csv("orphanet_all_parsed.txt", sep="\t", na_rep="NA", index=False,  encoding='utf-8')

