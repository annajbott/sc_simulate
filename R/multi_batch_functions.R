# Functions needed for splatter namespace

###########################################################################################

# Rewritten batch effects
splatSimBatchEffects_multi <- function(sim, params,
                                      nPatients = nPatients, patient.facLoc = patient.facLoc, patient.facScale = patient.facScale,
                                      nChannels = nChannels, channel.facLoc = channel.facLoc, channel.facScale = channel.facScale) {
    
    print("batch_batch_effects_multi")
    print("in")
    
    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    batch.facLoc <- getParam(params, "batch.facLoc")
    batch.facScale <- getParam(params, "batch.facScale")
    means.gene <- rowData(sim)$GeneMean
    
    for (idx in seq_len(nBatches)) {
        batch.facs <- getLNormFactors(nGenes, 1, 0.5, batch.facLoc[idx],
                                        batch.facScale[idx])
        batch.means.gene <- means.gene * batch.facs
        rowData(sim)[[paste0("BatchFacBatch", idx)]] <- batch.facs
    }
    
    patient.facLoc <- rep(patient.facLoc, nPatients*nBatches)
    patient.facScale <- rep(patient.facScale, nPatients*nBatches)
    

    for (idx in seq_len(nPatients*nBatches)) {
        patient.facs <- getLNormFactors(nGenes, 1, 0.5, patient.facLoc[idx],
                                        patient.facScale[idx])
        patient.means.gene <- means.gene * patient.facs
        
        rowData(sim)[[paste0("PatientFacPatient", idx)]] <- patient.facs

    }

    channel.facLoc <- rep(channel.facLoc, nChannels)
    channel.facScale <- rep(channel.facScale, nChannels)
    means.gene <- rowData(sim)$GeneMean

    for (idx in seq_len(nChannels)) {
        channel.facs <- getLNormFactors(nGenes, 1, 0.5, channel.facLoc[idx],
                                        channel.facScale[idx])
        channel.means.gene <- means.gene * channel.facs
        rowData(sim)[[paste0("ChannelFacChannel", idx)]] <- channel.facs
    }
    
    return(sim)
}

###########################################################################################

# Rewritten cell means

splatSimBatchCellMeans_multi <- function(sim, params) {

    print("batch_cell_means_multi")
    print("in")
    nBatches <- getParam(params, "nBatches")
    cell.names <- colData(sim)$Cell
    gene.names <- rowData(sim)$Gene
    gene.means <- rowData(sim)$GeneMean

    if (nBatches > 1) {
        batches <- colData(sim)$Batch
        batch.names <- unique(batches)

        batch.facs.gene <- rowData(sim)[, paste0("BatchFac", batch.names)]

        batch.facs.cell <- as.matrix(batch.facs.gene[,
                                                  as.numeric(factor(batches))])

        patient <- colData(sim)$Patient
        patient.names <- unique(patient)

        ########## here

        patient.facs.gene <- rowData(sim)[, paste0("PatientFac", patient.names)]
        patient.facs.cell <- as.matrix(patient.facs.gene[,
                                                  as.numeric(factor(patient))])

        channel <- colData(sim)$Channel
        channel.names <- unique(channel)

        channel.facs.gene <- rowData(sim)[, paste0("ChannelFac", channel.names)]
        channel.facs.cell <- as.matrix(channel.facs.gene[,
                                                  as.numeric(factor(channel))])

    } else {
        nCells <- getParam(params, "nCells")
        nGenes <- getParam(params, "nGenes")

        batch.facs.cell <- matrix(1, ncol = nCells, nrow = nGenes)
    }
    # Multiplied by all contributing to batch effect
    batch.means.cell <- patient.facs.cell * batch.facs.cell * channel.facs.cell * gene.means 

    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- gene.names

    # Captures all
    assays(sim)$BatchCellMeans <- batch.means.cell

    return(sim)
}

###########################################################################################
# Rewritten general splatSimulate

splatSimulate_multi_batches <- function(params = newSplatParams(),
                          method = c("single", "groups", "paths"),
                          verbose = TRUE, nChannels = 70,  nPatients = 14, 
                                        channel.facLoc = 0.05, channel.facScale = 0.05, 
                                        patient.facLoc = 0.15, patient.facScale = 0.15, ...) {
    
    print("simulate_multi_batches")
    # nPatients (number of patients per pool)
    checkmate::assertClass(params, "SplatParams")

    method <- match.arg(method)

    if (verbose) {message("Getting parameters...")}
    params <- setParams(params, ...)
    params <- expandParams(params)
    validObject(params)

    # Set random seed
    seed <- getParam(params, "seed")
    set.seed(seed)

    # Get the parameters we are going to use
    nCells <- getParam(params, "nCells")
    nGenes <- getParam(params, "nGenes")
    nBatches <- getParam(params, "nBatches")
    # Used here for pool
    batch.cells <- getParam(params, "batchCells")
    nGroups <- getParam(params, "nGroups")
    group.prob <- getParam(params, "group.prob")

    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to single mode")
        method <- "single"
    }

    if (verbose) {message("Creating simulation object...")}
    # Set up name vectors
    cell.names <- paste0("Cell", seq_len(nCells))
    gene.names <- paste0("Gene", seq_len(nGenes))
    batch.names <- paste0("Batch", seq_len(nBatches))
    
    # Names for my extra batches
    patient.names <- paste0("Patient", seq_len(nPatients*nBatches))
    channel.names <- paste0("Channel", seq_len(nChannels))



    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else if (method == "paths") {
        group.names <- paste0("Path", seq_len(nGroups))
    }

    # Create SingleCellExperiment to store simulation
    cells <-  data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    features <- data.frame(Gene = gene.names)
    rownames(features) <- gene.names
    sim <- SingleCellExperiment(rowData = features, colData = cells,
                                metadata = list(Params = params))

    # Make batches vector which is the index of param$batchCells repeated
    # params$batchCells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {rep(i, b[i])},
                      b = batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]
    
    # Repeat channel number, e.g. 140 channel 1, 140 channel 2, 140 channel 3... etc. 
    channels = rep((1:nChannels), each  = (nCells / nChannels))
    
    # Number of cells in each pool
    length_pool = nCells / nBatches
    
    # Patient label. Set of patients in each pool, e.g. 14 patients in pool1, patients 15-28 in pool 2. Randomised for each cell in that bracket. 
    patients = sapply(1:10, function(x) sample((((x*nPatients)-(nPatients-1)):(x*nPatients)), length_pool ,replace = TRUE)) %>% as.vector()

    # Channel and patient coldata 
    colData(sim)$Patient <- patient.names[patients]
    colData(sim)$Channel <- channel.names[channels]

    if (method != "single") {
        groups <- sample(seq_len(nGroups), nCells, prob = group.prob,
                         replace = TRUE)
        colData(sim)$Group <- factor(group.names[groups], levels = group.names)
    }

    if (verbose) {message("Simulating library sizes...")}
    sim <- splatSimLibSizes(sim, params)
    if (verbose) {message("Simulating gene means...")}
    sim <- splatSimGeneMeans(sim, params)
    if (nBatches > 1) {
        if (verbose) {message("Simulating batch effects...")}
        # Added extra parameter options
        sim <- splatSimBatchEffects(sim, params,
                                   nPatients = nPatients, patient.facLoc = patient.facLoc, patient.facScale = patient.facScale,
                                      nChannels = nChannels, channel.facLoc = channel.facLoc, channel.facScale = channel.facScale)
    }
    sim <- splatSimBatchCellMeans(sim, params)
    if (method == "single") {
        sim <- splatSimSingleCellMeans(sim, params)
    } else if (method == "groups") {
        if (verbose) {message("Simulating group DE...")}
        sim <- splatSimGroupDE(sim, params)
        if (verbose) {message("Simulating cell means...")}
        sim <- splatSimGroupCellMeans(sim, params)
    } else {
        if (verbose) {message("Simulating path endpoints...")}
        sim <- splatSimPathDE(sim, params)
        if (verbose) {message("Simulating path steps...")}
        sim <- splatSimPathCellMeans(sim, params)
    }
    if (verbose) {message("Simulating BCV...")}
    sim <- splatSimBCVMeans(sim, params)
    if (verbose) {message("Simulating counts...")}
    sim <- splatSimTrueCounts(sim, params)
    if (verbose) {message("Simulating dropout (if needed)...")}
    sim <- splatSimDropout(sim, params)

    if (verbose) {message("Done!")}
    return(sim)
}

###########################################################################################

