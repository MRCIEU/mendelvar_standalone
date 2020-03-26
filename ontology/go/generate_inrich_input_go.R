library("ontologyIndex")
library(tools)
library("hash")

args = commandArgs(trailingOnly=TRUE)

go <- args[1]
gaf <- args[2]
my_output <- args[3]

get_relation_names(go)
ontology <- get_ontology(go, propagate_relationships=c("is_a", "part_of"))
check(ontology, stop_if_invalid = FALSE)

ontology_slim <- get_ontology("goslim_generic.obo", propagate_relationships=c("is_a", "part_of"))
ontology_slim_ids <- ontology_slim$id
names(ontology_slim_ids) <- NULL
ontology_slim_ids <- gsub("GO:", "", ontology_slim_ids)
write(ontology_slim_ids, file = "goslim_generic.ids",
      ncolumns = 1,
      append = FALSE, sep = "\n")

my_data <- read.delim(gaf, stringsAsFactors = F, header=F, sep="\t",comment.char = '!', colClasses="character")

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
counter = 0
for (row in 1:nrow(my_data)) {
  counter = counter + 1
  print(counter)
  my_gene_id <- my_data[row,2]
  my_hpo <- my_data[row,5][1]
  if (has.key(my_hpo, h)) {
    my_ancestors <- values(h, keys=my_hpo)
    for (a in my_ancestors) {
    my_set_name <- ontology$name[a]
    names(my_set_name) <- NULL
    a <- gsub("GO:", "", a)
    my_id <- paste("gene:", my_gene_id, "GO:", a)
    g[my_id] <- data.frame(my_gene_id, a, my_set_name[1], stringsAsFactors = FALSE)
    }
  }}

f <- do.call(rbind, values(g))
gene <- f[seq(1, length(f), 3)]
hpo <- f[seq(2, length(f), 3)]
term <- f[seq(3, length(f), 3)]
final_df <- data.frame("gene_id"=gene, "gene_set"=hpo, "name"=term, stringsAsFactors=FALSE)

#Remove duplicates from data frame
#Remove 0005575 (cellular_component), 0008150 (biological_process) and 0003674 (molecular_function)
final_df <- final_df[final_df$gene_set != "0005575" & final_df$gene_set != "0008150" & final_df$gene_set != "0003674",]
final_df <- unique(final_df)
#Write to file
write.table(final_df, my_output, sep="\t", col.names=F, row.names=F, quote=F)

#Filter to keep only GO slim terms.
slim_out <- final_df[final_df$gene_set %in% ontology_slim_ids,]
write.table(slim_out, "go_slim_inrich.txt", sep="\t", col.names=F, row.names=F, quote=F)