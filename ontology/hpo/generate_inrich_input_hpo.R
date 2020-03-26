library("ontologyIndex")
library("hash")
library("tools")
args = commandArgs(trailingOnly=TRUE)

hpo <- args[1]
gaf <- args[2]
my_output <- args[3]

get_relation_names(hpo)
ontology <- get_ontology(hpo, propagate_relationships=c("is_a", "part_of"))
check(ontology, stop_if_invalid = FALSE)

my_data <- read.delim(gaf, stringsAsFactors = F, header=F, sep="\t",comment.char = '!', colClasses="character")

#Add HP
my_data <- transform(my_data, V5=sprintf('HP:%s', V5))

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
    a <- gsub("HP:", "", a)
    my_id <- paste("gene:", my_gene_id, "_HP:", a)
    g[my_id] <- data.frame(my_gene_id, a, my_set_name[1], stringsAsFactors = FALSE)
    }
  }}

f <- do.call(rbind, values(g))
gene <- f[seq(1, length(f), 3)]
hpo <- f[seq(2, length(f), 3)]
term <- f[seq(3, length(f), 3)]
final_df <- data.frame("gene_id"=gene, "gene_set"=hpo, "name"=term, stringsAsFactors=FALSE)


#Remove duplicates from data frame
#Remove 0000001 (All) and 0000118 (phenotypic abnormality)
final_df <- final_df[final_df$gene_set != "0000001" & final_df$gene_set != "0000118",]
final_df <- unique(final_df)
#Write to file
write.table(final_df, my_output, sep="\t", col.names=F, row.names=F, quote=F)
