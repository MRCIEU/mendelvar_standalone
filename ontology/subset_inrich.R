library(tools)

args = commandArgs(trailingOnly=TRUE)

g_d <- args[1]
inrich <- args[2]
inrich_out_file <- args[3]

g_d_df <- read.delim(g_d, stringsAsFactors = F, header=T, sep="\t", colClasses="character")
inrich_df <- read.delim(inrich, stringsAsFactors = F, header=F, sep="\t", colClasses="character")
names(inrich_df) <- c("gene_id", "gene_set", "name")
inrich_out <- inrich_df[inrich_df$gene_id %in% g_d_df$hgnc_id,]
write.table(inrich_out, inrich_out_file, sep="\t", col.names=F, row.names=F, quote=F)