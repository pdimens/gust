suppressPackageStartupMessages(library(ggtree))
pdf(NULL)
input <- read.tree(snakemake@input[[1]])
ggtree(input) + geom_tiplab() + geom_nodelab(size=3, col="dodgerblue")
ggsave(snakemake@output[[1]])
