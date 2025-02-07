library(plyr)
library(tidyverse)
library(splatter)
library(Matrix)
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

nCells = 700000
nChannels = 70
nBatches = 10
nPatients = 140 # 14 per pool, 140 in total
nGenes = 4000
distinct_barcode_pool <- FALSE # TRUE means for each pool, randomly sample barcodes. If FALSE distinct barcodes across each cell and pool (if nCells < 737,280)

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


## Read in highest expressed PBMC ensembl gene names 
## Replace Gene1, Gene2 ... with Ensembl names

highest_expressed_list <-  (sce_pbmc %>% counts) %>% rowSums() %>% order() %>% rev() # Highest expressed genes from PBMCs
highest_expressed_list_top <- highest_expressed_list[1:nGenes]
highest_expressed_genes <- sce_pbmc[highest_expressed_list_top,] %>% rownames
rownames(sim) <- highest_expressed_genes




#########

### Read barcodes in
barcodes <- read_delim("../data/737K-august-2016.txt", delim = "\t", col_names =  FALSE)
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

if(!distinct_barcode_pool & nCells < (barcodes %>% unique %>% length)){
cell_barcodes_all <- sample(barcodes, nCells, replace = FALSE)}

for(idx in seq(nBatches)){

    # Subset large matrix into each pool
    start = (idx*cells_in_pool) - (cells_in_pool - 1) # Start of each pool
    end = idx*cells_in_pool # End of each pool

    sce_subset = sim[,start:end]
    metadata[[idx]] <- as.data.frame(sce_subset@colData)
    
    # Distinct CBs for whole experiment or for each pool only
    if(!distinct_barcode_pool & nCells < (barcodes %>% unique %>% length)){
    cell_barcodes <- cell_barcodes_all[start:end]
    } else{
    cell_barcodes <- sample(barcodes, cells_in_pool, replace = FALSE)
    }

    colData(sce_subset)$Cell <- cell_barcodes
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
        barcodes_file <- paste0(c(dir, "/quants_mat_cols.txt"), collapse ="" )
        gene_file <- paste0(c(dir, "/quants_mat_rows.txt"), collapse ="" )
        count_file <- paste0(c(dir, "/quants_mat.csv"), collapse ="" )

        # Write to files
        write.table(rownames(sce_channel_subset), file= gene_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(colnames(sce_channel_subset), file= barcodes_file, quote=FALSE, col.names=FALSE, row.names=FALSE)
        write.table(counts(sce_channel_subset), file= count_file, quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",") 

    }
}


# Could do with if statements in for loop, but would slow down considerably
metadata = do.call(rbind, metadata)


#test <- read.csv("out_dir/channel_1/quants_mat.csv")
#example_mat <- read.csv("out_dir/channel_1/quants_mat.csv", header = FALSE)


mat_all <- list() 
cols_all <- list()
rows_all <- list()

for(idx in seq(nChannels)){
    print(idx)
    file_mat <- paste0("out_dir/channel_", idx, "/quants_mat.csv")
    file_cols <- paste0("out_dir/channel_", idx, "/quants_mat_cols.txt")
    file_rows <- paste0("out_dir/channel_", idx, "/quants_mat_rows.txt")
    
    mat <- read.csv(file = file_mat, header = FALSE)
    cols <- read.csv(file = file_cols, header = FALSE)
    rows <- read.csv(file = file_rows, header = FALSE)
    
    
    cols <- cols$V1
    rows <- rows$V1
    
    colnames(mat) <- cols
    rownames(mat) <- rows
    
    mat_all[[idx]] <- mat
    cols_all[[idx]] <- cols
    rows_all[[idx]] <- rows
    
}

mat_all = do.call(cbind, mat_all)

rownames(mat_all) <- rows
cols_all = do.call(rbind, cols_all)
rows_all = do.call(rbind, rows_all)

# Replace names in group with sample name
metadata$Group <- revalue(metadata$Group, c("Group1"="Severe_COVID", "Group2"="HospMild_COVID", "Group3"="Flu", "Group4"="Healthy", "Group5"="HcwMild_COVID"))



mat_name <- paste0(nCells, "_", "ch", nChannels, "_", "final_mat.csv")
col_name <- paste0(nCells, "_", "ch", nChannels, "_", "final_cols.txt")
row_name <- paste0(nCells, "_", "ch", nChannels, "_", "final_rows.txt")
write.csv(mat_all, file = mat_name)

write.table(cols_all, file=col_name, row.names = FALSE, col.names = FALSE)
write.table(rows_all, file=row_name, row.names = FALSE, col.names = FALSE)

meta_name <- paste0(nCells, "_", "ch", nChannels, "_", "metadata.csv")
write.csv(metadata, file = meta_name)

# Generate market matrix

sparse.all <- Matrix(as.matrix(mat_all), sparse = T)

mat_name <- paste0(nCells, "_", "ch", nChannels, "_", "matrix.mtx")
col_name <- paste0(nCells, "_", "ch", nChannels, "_", "barcodes.tsv")
row_name <- paste0(nCells, "_", "ch", nChannels, "_", "genes.tsv")
writeMM(obj=sparse.all, file = mat_name)

write.table(colnames(mat_all), file=col_name, row.names = FALSE, col.names = FALSE)
write.table(rownames(mat_all), file=row_name, row.names = FALSE, col.names = FALSE)
