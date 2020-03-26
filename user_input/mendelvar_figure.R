library(stringr)
library(ggplot2)
library(scales)
library(viridis)
library(tools)

args = commandArgs(trailingOnly=TRUE)
input_file <- args[1]
inrich_parsed <- read.delim(input_file, stringsAsFactors = F, header=T, sep="\t")
if (nrow(inrich_parsed) >= 1)
{
inrich_parsed$number_genes_overlap <- str_count(inrich_parsed$Overlapping_genes_hgnc_id, ",") + 1
inrich_parsed$gene_ratio_term <- inrich_parsed$number_genes_overlap / inrich_parsed$Gene_set_size 
inrich_parsed$log_pval <- -log10(inrich_parsed$Empirical_p_value)
inrich_parsed_subset <- inrich_parsed[1:50,]
# lock in factor level order
#Need to shorten labels on Gene set description
Encoding(inrich_parsed_subset$Gene_set_description) <- "latin1"
inrich_parsed_subset$Gene_set_description <- substr(inrich_parsed_subset$Gene_set_description, 1, 68)
inrich_parsed_subset$Gene_set_description <- factor(inrich_parsed_subset$Gene_set_description, levels = unique(inrich_parsed_subset$Gene_set_description))

basic <- ggplot(inrich_parsed_subset) + geom_point(aes(x = inrich_parsed_subset$number_genes_overlap, y = inrich_parsed_subset$Gene_set_description, size = inrich_parsed_subset$gene_ratio_term, color = inrich_parsed_subset$Empirical_p_value)) + scale_color_viridis(option = "D", direction=-1) + scale_y_discrete(limits = rev(levels(inrich_parsed_subset$Gene_set_description)))

final <- basic + expand_limits(x = 0, y = 0) + scale_x_continuous(breaks=pretty_breaks()) + labs(color = "Empirical p-value", size="Ratio of gene overlap") + theme(panel.grid.major.y = element_line(colour = "gray",size=0.2,  linetype="dotted"), panel.grid.major.x = element_line(NA), panel.grid.minor = element_line(NA), panel.background = element_rect(fill=NA)) + xlab("Number of genes in the overlap") + ylab("") + theme(axis.line.x = element_line(color="gray", size = 0.1), legend.key=element_blank())
} else {
#Empty plot if no output
final <- ggplot() + theme_void()
}
outfile_prefix <- strsplit(input_file, ".out.inrich.parsed")[[1]]

my_pdf <- paste(outfile_prefix, ".pdf", sep="")
my_png <- paste(outfile_prefix, ".png", sep="")

ggsave(my_pdf, final, dpi=300, height=7.21, width=8.32)