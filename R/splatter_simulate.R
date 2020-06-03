library(tidyverse)
library(splatter)
library(scater)
library(optparse)

option_list <- list(make_option(c("--sce"), default="must_specify",
	     help="The sce rds file to estimate parameters from"),
	     make_option(c("--ngenes"), default="all",
	     help="Number of genes to simulate, default does not subset"))

opt <- parse_args(OptionParser(option_list=option_list))
#sce <- read.table(opt$sce,header = TRUE, fill = T)
subset_genes <- opt$ngenes
#read.table(opt$ngenes,header = TRUE, fill=T)

# Load SCE from 10X PBMC data (1000 cells)- processed with salmon alevin
sce_pbmc <- readRDS("ground_truth_1000_sce.rds")

# Estimate paramters for splatter
params_pmbc <- splatEstimate(sce_pbmc)

# Set number of cells to 8000s
ncells <- 400000
params_pmbc <- setParam(params_pmbc, "batchCells", ncells)

# Set number of genes (if we want a subset)
if(subset_genes == "all"){
ngenes <- length(rownames(sim))
}
ngenes <- opt$ngenes
params_pmbc <- setParam(params_pmbc, "nGenes", ngenes)

nMat <-10

dir.create(path = "out_single_platter_dir")

for(i in 1:nMat){
# Simulated SCE
# Change seed for each simulation, otherwise deterministic
params_pmbc <- setParam(params_pmbc, "seed", i*20)
sim <- splatSimulate(params_pmbc)
# Get counts with counts(sim)

# Currently just gene1,2,3 etc and cell1,2,3 etc. Can give ensembl gene names
# gene_names <- rownames(sce_pbmc)[1:ngenes]

dir <- paste0(c("out_single_platter_dir/", "channel_", i), collapse = "")
dir.create(path = dir)

barcodes_file <- paste0(c(dir, "/quants_mat_rows.txt"), collapse ="" )
gene_file <- paste0(c(dir, "/quants_mat_cols.txt"), collapse ="" )
count_file <- paste0(c(dir, "/quants_mat.csv"), collapse ="" )

write.table(rownames(sim), file= gene_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(colnames(sim), file= barcodes_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= count_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",") 
rm(sim)

}
