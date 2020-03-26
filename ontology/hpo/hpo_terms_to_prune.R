library("ontologyIndex")
library(tools)

args = commandArgs(trailingOnly=TRUE)

hpo <- args[1]
disease <- args[2]
my_output <- args[3]

get_relation_names(hpo)
ontology <- get_ontology(hpo, propagate_relationships=c("is_a", "part_of"))
check(ontology, stop_if_invalid = FALSE)
to_remove <- get_descendants(ontology, roots=c("HP:0000005", "HP:0031797", "HP:0040279", "HP:0012823"))
to_remove <- gsub("HP:", "", to_remove)
my_data <- read.delim(disease, stringsAsFactors = F, header=T, sep="\t", colClasses="character")
hpos <- my_data$hpo
a <- strsplit(hpos, ";", fixed = TRUE, perl = FALSE, useBytes = FALSE)
#Test if removal in the set really works. Also tested on raw file.
a[[1]] <- c(a[[1]], "0000005")
new <- lapply(1:length(a), function(n) setdiff(a[[n]], to_remove))
new2 <- lapply(1:length(new), function(n) paste(new[[n]], collapse=";"))
my_data$hpo <- as.character(new2)
write.table(my_data, my_output, sep="\t", col.names=T, row.names=F, quote=F)