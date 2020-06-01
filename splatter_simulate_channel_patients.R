library(tidyverse)
library(splatter)
#library(optparse)
source('multi_batch_functions.R')


# Rewritten certain functions in Splat simulate 
# Loaded from splatter_methods.R
# Import them into namespace

#######################
# Rename in namespace #
#######################

# Batch effect
environment(splatSimBatchEffects_multi) <- asNamespace('splatter')
assignInNamespace("splatSimBatchEffects", splatSimBatchEffects_multi, ns = "splatter")

# Batch cell means
environment(splatSimBatchCellMeans_multi) <- asNamespace('splatter')
assignInNamespace("splatSimBatchCellMeans", splatSimBatchCellMeans_multi, ns = "splatter")

# Splat simulate multi batches
environment(splatSimulate_multi_batches) <- asNamespace('splatter')
assignInNamespace("splatSimulate", splatSimulate_multi_batches, ns = "splatter")


### Maybe optparse feature for nGenes etc.

nCells = 1400
nChannels = 10
nBatches = 10
nPatients = 14 # 14 per pool, 140 in total
nGenes = 8000

# Number of cells in each batch (pool)
batchcells = rep(nCells/nBatches, nBatches)

# Load SCE from 10X PBMC data (1000 cells)- processed with salmon alevin
sce_pbmc <- readRDS("ground_truth_1000_sce.rds")
message("SCE loaded")

# Estimate paramters for splatter
params_pmbc <- splatEstimate(sce_pbmc)
message("Parameters estimated")

# Set number of cells to 1.4 million (round number with 70 channels) - not needed
#ncells <- 1400000
#params_pmbc <- setParam(params_pmbc, "batchCells", ncells)
params_pmbc <- setParam(params_pmbc, "lib.loc", 11)
params_pmbc <- setParam(params_pmbc, "nGenes", nGenes)

# Set seed
set.seed(1234)


# 140,000 cells in each batch. Equally likely to be in one of 5 groups (random)
# Batch = pool
# 70 channels
# 14 patients (140 total, 14 in each pool)
# facLoc and facScale affects mean and sd of batch/group effects 
sim <- splatSimulate_multi_batches(params = params_pmbc,
                            batchCells = batchcells, method = "groups", group.prob = c(0.2,0.2,0.2,0.2,0.2),
                            verbose = TRUE, nChannels = nChannels,  nPatients = nPatients, 
                            channel.facLoc = 0.05, channel.facScale = 0.05, 
                            patient.facLoc = 0.15, patient.facScale = 0.15, 
                            de.prob = c(0.1, 0.1, 0.1, 0.2, 0.2), de.facLoc = 0.2, de.facScale = 0.4,)

## Read in random ensembl gene names here
## Replace Gene1, Gene2 ... with Ensembl names
# .. 
# ...
# ....

#########

### Read barcodes in
barcodes <- read_delim("737K-august-2016.txt", delim = "\t", col_names =  FALSE)
barcodes <- barcodes$X1


#  There are 737280 barcodes < number of cells -- could randomly sample from each pool- as won't be in the same run anyway
# cell_barcodes <- sample(barcodes,ncells, replace = FALSE )

# Split up cells into pools, then give barcodes, then split into channels?

dir.create(path = "out_dir")

cells_in_pool = nCells/nBatches
channels_per_pool = nChannels/nBatches
cells_in_channel = cells_in_pool/channels_per_pool
channel_no = 0 

metadata <-  list() 

for(idx in seq(nBatches)){

    # Subset large matrix into each pool
    start = (idx*cells_in_pool) - (cells_in_pool - 1) # Start of each pool
    end = idx*cells_in_pool # End of each pool

    sce_subset = sim[,start:end]
    metadata[[idx]] <- as.data.frame(sce_subset@colData)
    
    cell_barcodes <- sample(barcodes, cells_in_pool, replace = FALSE)

    colData(sce_subset)$Cell <- cell_barcodes
    # rownames too
    rownames(colData(sce_subset)) <- cell_barcodes

    # Each channel per pool
    for(chnl in seq(channels_per_pool)){
        channel_no = channel_no + 1
        start = (chnl*cells_in_channel) - (cells_in_channel - 1)
        end = chnl*cells_in_channel # End of each pool
        sce_channel_subset = sce_subset[,start:end]
        # Create directories
        dir <- paste0(c("out_dir/", "channel_", channel_no), collapse = "")
        dir.create(path = dir)
        barcodes_file <- paste0(c(dir, "/quants_mat_rows.txt"), collapse ="" )
        gene_file <- paste0(c(dir, "/quants_mat_cols.txt"), collapse ="" )
        count_file <- paste0(c(dir, "/quants_mat.csv"), collapse ="" )

        # Write to files
        write.table(rownames(sce_channel_subset), file= gene_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(colnames(sce_channel_subset), file= barcodes_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(counts(sce_channel_subset), file= count_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",") 

    }
}

metadata = do.call(rbind, metadata)


#test <- read.csv("out_dir/channel_1/quants_mat.csv")
example_mat <- read.csv("out_dir/channel_1/quants_mat.csv", header = FALSE)


mat_all <- list() 


for(idx in seq(nBatches)){
    print(idx)
    file_mat <- paste0("out_dir/channel_", idx, "/quants_mat.csv")
    file_cols <- paste0("out_dir/channel_", idx, "/quants_mat_cols.txt")
    file_rows <- paste0("out_dir/channel_", idx, "/quants_mat_rows.txt")
    
    mat <- read.csv(file = file_mat, header = FALSE)
    cols <- read.csv(file = file_cols, header = FALSE)
    rows <- read.csv(file = file_rows, header = FALSE)
    
    
    cols <- cols$V1
    rows <- rows$V1
    
    colnames(mat) <- rows
    rownames(mat) <- cols
    
    mat_all[[idx]] <- mat
    
}

mat_all = do.call(cbind, mat_all)

write.csv(mat_all, file = "final_mat.csv")
