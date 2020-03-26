library("hash")
library("ontologyIndex")
library(tools)

args = commandArgs(trailingOnly=TRUE)

do <- args[1]
gaf <- args[2]
my_output <- args[3]

get_relation_names(do)
ontology <- get_ontology(do, propagate_relationships=c("is_a", "part_of"))
check(ontology, stop_if_invalid = FALSE)

#Get ids of terms in ontology slim for later filtering.
ontology_slim <- get_ontology("DO_AGR_slim.obo", propagate_relationships=c("is_a", "part_of"))
ontology_slim_ids <- ontology_slim$id
names(ontology_slim_ids) <- NULL
ontology_slim_ids <- gsub("DOID:", "", ontology_slim_ids)
write(ontology_slim_ids, file = "DO_AGR_slim.ids",
      ncolumns = 1,
      append = FALSE, sep = "\n")

my_data <- read.delim(gaf, stringsAsFactors = F, header=F, sep="\t",comment.char = '!', colClasses="character")

#Add HP
my_data <- transform(my_data, V5=sprintf('DOID:%s', V5))

#Prepare a file ready to be used by INRICH:  gene ID (using hgnc ID), gene-set ID, gene-set name
h <- hash()
g <- hash()
#Save all the hpo - ancestors relationships to a hash. Ancestors include the query hpo terms itself too.
for (row in 1:nrow(my_data)) {
  my_hpo <- my_data[row,5][1]
  my_ancestors <- get_ancestors(ontology, my_hpo)
  if (has.key(my_hpo, h)) {
  }
  else{
    h[my_hpo] <- my_ancestors
  }
  
  }

#Save input for INRICH: gene name, gene set id, gene set name
for (row in 1:nrow(my_data)) {
  my_gene_id <- my_data[row,2]
  my_hpo <- my_data[row,5][1]
  if (has.key(my_hpo, h)) {
    my_ancestors <- values(h, keys=my_hpo)
    for (a in my_ancestors) {
    my_set_name <- ontology$name[a]
    names(my_set_name) <- NULL
    a <- gsub("DOID:", "", a)
    my_id <- paste("gene:", my_gene_id, "_DOID:", a)
    g[my_id] <- data.frame(my_gene_id, a, my_set_name[1], stringsAsFactors = FALSE)
    }
  }}

f <- do.call(rbind, values(g))
gene <- f[seq(1, length(f), 3)]
hpo <- f[seq(2, length(f), 3)]
term <- f[seq(3, length(f), 3)]
final_df <- data.frame("gene_id"=gene, "gene_set"=hpo, "name"=term, stringsAsFactors=FALSE)

#Remove duplicates from data frame
#Remove 7 (disease), 630 (genetic disease) and 7 (disease of anatomical entity)
final_df <- final_df[final_df$gene_set != "4" & final_df$gene_set != "7" & final_df$gene_set != "630",]
final_df <- unique(final_df)
#Write to file
write.table(final_df, my_output, sep="\t", col.names=F, row.names=F, quote=F)

#Filter to keep only DO slim terms.
slim_out <- final_df[final_df$gene_set %in% ontology_slim_ids,]
write.table(slim_out, "do_slim_inrich.txt", sep="\t", col.names=F, row.names=F, quote=F)