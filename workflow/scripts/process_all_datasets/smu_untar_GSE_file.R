# Script to process the raw dataset from the Sean M Wu Paper

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/smu_untar_GSE_file.rda")
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Opening log file to collect all messages, warnings and errors
message("Opening log file")
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

### LOADING FILES ============================================================

# Load in the necessary packages
message("Loading packages")
suppressPackageStartupMessages({
  library(hdf5r)
  library(Seurat)
  library(DropletUtils)
  library(SoupX)
  library(tidyverse)
  library(SeuratWrappers)
})


### UNPACKING THE TAR FILE ===================================================

message("Unpacking the tar file")

# Dynamically setting exdir based on tar file location
tar_file <- snakemake@input$tar_file
output_dir <- dirname(tar_file)

# Unpack the tar file into the dynamically created directory
untar(tar_file, exdir = output_dir)
message("Tar file unpacked successfully")

# Get the list of files in the directory for further processing
untarred_files <- list.files(output_dir, pattern = "\\.h5$", full.names = TRUE)
message("H5 files extracted: ", paste(basename(untarred_files), collapse = ", "))

### LOAD IN THE OBJECTS =======================================================

# List of raw and filtered files
raw_files <- untarred_files[grepl("raw_feature_bc_matrix.h5$", untarred_files)]
filtered_files <- untarred_files[grepl("filtered_feature_bc_matrix.h5$", untarred_files)]

raw_data_list <- lapply(raw_files, Read10X_h5)
names(raw_data_list) <- basename(raw_files)

filtered_data_list <- lapply(filtered_files, Read10X_h5)
names(filtered_data_list) <- basename(filtered_files)


### ABMIENT RNA CORRECTION ====================================================

# Create SoupChannel objects
soup_channel_list <- lapply(names(raw_data_list), function(name) {
  filtered_name <- gsub("raw", "filtered", name)
  raw_data <- raw_data_list[[name]]$`Gene Expression`
  filtered_data <- filtered_data_list[[filtered_name]]$`Gene Expression`
  sc <- SoupChannel(tod = raw_data, toc = filtered_data)
  sc$sampleName <- name
  return(sc)
})
names(soup_channel_list) <- names(raw_data_list)

# Run preliminary clustering so that SoupX can calculate ambient RNA
preliminary_clustering <- function(sc) {
  # Create and cluster a Seurat object 
  seurat_obj <- CreateSeuratObject(sc$toc)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(object = seurat_obj, nfeatures = 2000, verbose = FALSE, selection.method = 'vst')
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  
  # Add the clusters to the SoupX object
  sc <- setClusters(sc, setNames(seurat_obj$seurat_clusters, colnames(sc$toc)))
  return(sc)
  
}
soups_grouped <- lapply(soup_channel_list, preliminary_clustering)

# Correct for Ambient RNA
run_soupx <- function(sc) {
  # Define output directory based on the combined_seurat output
  output_dir <- dirname(snakemake@output$processed_data)
  
  # Construct a filename for the diagnostic plot using the sample name
  diag_file <- file.path(output_dir, paste0(sc$sampleName, "_autoEstCont_diag.pdf"))
  
  # Open a PDF device to capture the diagnostic plot
  pdf(file = diag_file)
  sc <- autoEstCont(sc, doPlot = TRUE)  # This will now plot into the PDF device
  dev.off()
  
  # Now adjust counts using the corrected SoupChannel object
  adjusted_counts <- adjustCounts(sc, roundToInt = TRUE)
  return(adjusted_counts)
}
adjusted_counts_list <- lapply(soups_grouped, run_soupx)
names(adjusted_counts_list) <- names(soups_grouped)


### CREATE SEURAT OBJECTS WITH ADJUSTED COUNTS AND CMO MATRIX =================

# Create Seurat objects and include the CMO matrix
seurat_objects <- lapply(names(adjusted_counts_list), function(sample_name) {
  # Grab the sample data 
  adjusted_counts <- adjusted_counts_list[[sample_name]]
  filtered_data <- filtered_data_list[[gsub("raw", "filtered", sample_name)]]
  cmo_data <- filtered_data$`Antibody Capture`
  
  # Create Seurat object with adjusted counts and CMO data
  seurat_obj <- CreateSeuratObject(counts = adjusted_counts)
  seurat_obj[["HTO"]] <- CreateAssayObject(counts = cmo_data)
  return(seurat_obj)
})
names(seurat_objects) <- names(adjusted_counts_list)

### DEMULTIPLEXING ===========================================================
#### For demultiplexing, I ran into some errors because some of the cell hashes
#### had such few instances -> I ended up performing demultiplexing, exluding
#### hashes with less than 1e6 rowSums. While this seems harsh, manually checking
#### each hash showed that if the rowSum wasn't >1e6, it was <50

# Demultiplex with HTODemux
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  # Remove CMO hashtags with low counts - else run into an error
  m <- seurat_obj[["HTO"]]@counts
  m <- m[which(rowSums(m) > 1e6),]
  seurat_obj[['HTOnew']] <- CreateAssayObject(counts = m)
  seurat_obj <- NormalizeData(seurat_obj, assay = "HTOnew", normalization.method = "LogNormalize")
  
  # Run HTODemux
  seurat_obj <- HTODemux(seurat_obj, assay = "HTOnew")
  return(seurat_obj)
})

# Remove doublet and negative cells, keeping only singlets
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  # Set the identity classes to the HTO classification and subset for singlets
  Idents(seurat_obj) <- "HTOnew_classification.global"
  return(subset(seurat_obj, idents = "Singlet"))
})


### QUALITY CONTROL ==========================================================

# Calculate the percentage of mitochondrial & ribosomal genes
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  # Percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Percentage of ribosomal genes
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  
  return(seurat_obj)
})

# List of QC metrics
qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")

# Function to perform QC filtering
perform_qc_filtering <- function(seurat_obj) {
  # Extract QC metrics data
  qc_data <- seurat_obj@meta.data[, qc_metrics]
  
  # Calculate median and MAD for each QC metric
  qc_medians <- apply(qc_data, 2, median, na.rm = TRUE)
  qc_mads <- apply(qc_data, 2, mad, na.rm = TRUE)
  
  # Calculate lower and upper cutoffs
  lower_cutoffs <- qc_medians - 3 * qc_mads
  upper_cutoffs <- qc_medians + 3 * qc_mads
  
  # Ensure lower cutoffs are not negative
  lower_cutoffs[lower_cutoffs < 0] <- 0
  
  # Identify cells that pass all QC thresholds
  cells_to_keep <- rownames(seurat_obj@meta.data)
  
  for (metric in qc_metrics) {
    lower <- lower_cutoffs[metric]
    upper <- upper_cutoffs[metric]
    values <- seurat_obj@meta.data[[metric]]
    
    # Identify cells within the thresholds for this metric
    keep_cells_metric <- cells_to_keep[values >= lower & values <= upper]
    
    # Update cells_to_keep by intersecting with keep_cells_metric
    cells_to_keep <- intersect(cells_to_keep, keep_cells_metric)
  }
  
  # Subset the Seurat object to keep only high-quality cells
  seurat_obj_filtered <- subset(seurat_obj, cells = cells_to_keep)
  
  # Optionally, you can print the number of cells before and after filtering
  message("Number of cells before QC filtering: ", ncol(seurat_obj))
  message("Number of cells after QC filtering: ", ncol(seurat_obj_filtered))
  
  return(seurat_obj_filtered)
}

# Apply the QC filtering to each Seurat object in the list
seurat_objects <- lapply(seurat_objects, perform_qc_filtering)


### MERGE AND NORMALIZE RUNS =================================================

# Merge the Seurat objects into one combined object
combined_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[-1],
  add.cell.ids = c("Run1", "Run2", "Run3"),
  project = "Combined_Seurat"
)

# Normalize the RNA assay counts
combined_seurat <- NormalizeData(combined_seurat, assay = "RNA")
combined_seurat <- JoinLayers(combined_seurat, overwrite = TRUE)
combined_seurat <- ScaleData(combined_seurat)


### IDENTIFY VARIABLE FEATURES AND INTEGRATE BATCHES =========================

# Add "run" and "individual" to metadata
combined_seurat$run <- sapply(strsplit(colnames(combined_seurat), "_"), `[`, 1)
combined_seurat$individual <- sapply(strsplit(combined_seurat$HTOnew_maxID, "-"), `[`, 1)
message("Added 'individual' and 'run' information to metadata")

# Split the Seurat object into a list by 'run' (dataset) and identify variable features in each run
message("Identifying variable features in each run")
seurat_list <- SplitObject(combined_seurat, split.by = "run")
seurat_list <- lapply(seurat_list, FindVariableFeatures)

# Select integration features across datasets
message("Selecting integration features")
features <- SelectIntegrationFeatures(object.list = seurat_list)

# Perform integration using RunFastMNN, specifying 'individual' as the batch variable
message("Running RunFastMNN")
combined_seurat <- RunFastMNN(object.list = seurat_list, 
                              features = features, 
                              batch = 'individual')
                             
message("Performed integration using RunFastMNN with 'individual' as batch variable")


### PERFORM CLUSTERING =======================================================

# After RunFastMNN, skip separate PCA and use MNN reduction consistently
message("Constructing the nearest neighbors graph")
combined_seurat <- FindNeighbors(combined_seurat, 
                                 reduction = "mnn",  # Use MNN reduction
                                 dims = 1:20)

message("Running UMAP for dimensionality reduction")
combined_seurat <- RunUMAP(combined_seurat, 
                           dims = 1:20,
                           reduction = "mnn",  # Use MNN reduction
                           reduction.name = "umap",
                           reduction.key = "UMAP_")

# Perform unsupervised clustering
message("Performing unsupervised clustering")
combined_seurat <- FindClusters(combined_seurat, resolution = 0.5)

# Visualize UMAP plot colored by clusters
message("Visualizing clusters on UMAP plot")
umap_clusters <- DimPlot(combined_seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +  NoLegend() + ggtitle("UMAP: Clusters")
umap_sample <- DimPlot(combined_seurat, reduction = "umap", group.by = "individual", label = TRUE) +  NoLegend() + ggtitle("UMAP: Individual IDs")

# Save UMAP plots
output_dir <- dirname(snakemake@output$processed_data)
ggsave(filename = file.path(output_dir, "umap_clusters.pdf"),
       plot = umap_clusters, device = "pdf", width = 6, height = 5)
ggsave(filename = file.path(output_dir, "umap_sample.pdf"),
       plot = umap_sample, device = "pdf", width = 6, height = 5)


### IDENTIFYING CARDIOMYOCYTE SPECIFIC CLUSTER ===============================

# Visualize key cardiomyocyte markers
cardiomyocyte_markers <- c("TNNT2", "MYL2", "MYH7", "TNNI1")
p1 <- FeaturePlot(combined_seurat, 
            features = cardiomyocyte_markers,
            ncol = 2,
            order = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90')

# Additional cardiac transcription factors/markers
cardiac_tfs <- c("NKX2-5", "IRX4", "TBX5", "HCN4")
p2 <- FeaturePlot(combined_seurat,
            features = cardiac_tfs,
            ncol = 2,
            order = TRUE,
            min.cutoff = 'q10',
            max.cutoff = 'q90')

# Specific markers that Liangfu sent
cardiac_tfs <- c("TNNT2", "MYH6", "MYH7", "ACTC1", "TPM1")
liangfu_markers <- FeaturePlot(combined_seurat, 
                               features = cardiac_tfs, 
                               ncol = 3, 
                               order = TRUE, 
                               min.cutoff = 'q10', 
                               max.cutoff = 'q90')

# Save the FeaturePlots to the same output directory
ggsave(filename = file.path(output_dir, "cardiomyocyte_markers.pdf"),
       plot = p1, device = "pdf", width = 8, height = 6)
ggsave(filename = file.path(output_dir, "additional_cardiomyocyte_tfs.pdf"),
       plot = p2, device = "pdf", width = 8, height = 6)
ggsave(filename = file.path(output_dir, "liangfu_markers.pdf"),
       plot = liangfu_markers, device = "pdf", width = 8, height = 6)


### COLORING UMAP BY SAMPLE ===================================================

# Create UMAP colored by hash.ID and assign it to a variable
hash_plot <- DimPlot(combined_seurat, 
                     reduction = "umap", 
                     group.by = "hash.ID", 
                     pt.size = 0.5, 
                     label = TRUE) + 
  ggtitle("Cell Distribution by Hash ID")

# Define output directory based on processed data output (if not already defined)
output_dir <- dirname(snakemake@output$processed_data)

# Save the hash.ID UMAP plot as PDF
ggsave(filename = file.path(output_dir, "umap_hashID.pdf"),
       plot = hash_plot, device = "pdf", width = 6, height = 5)


# Plot days 15 and 30 for each sample with a custom color vector
color_vector <- rep("lightgrey", length(unique(combined_seurat$hash.ID)))
names(color_vector) <- unique(combined_seurat$hash.ID)
color_vector["WTC-DAY30"] <- "#F8766D"
color_vector["WTC-DAY15"] <- "#7CAE00"
color_vector["LMNA-DAY15"] <- "#00BFC4"
color_vector["LMNA-DAY30"] <- "#C77CFF"

days_plot <- DimPlot(combined_seurat, 
                     reduction = "umap",
                     group.by = "hash.ID",
                     cols = color_vector,
                     pt.size = 0.5) +
  ggtitle("Days 15 and 30 Cells")

# Save the Days plot as PDF
ggsave(filename = file.path(output_dir, "umap_days15_30.pdf"),
       plot = days_plot, device = "pdf", width = 6, height = 5)


### SUBSETTING CARDIOMYOCYTES =================================================

# So it looks like clusters 12 and 14 are the Day 30 WTC and LMNA cells with highest expression of cardiomyocyte markers
# cool

# Subset clusters 12 and 14
cardiomyocyte_clusters <- subset(combined_seurat, seurat_clusters %in% c(12, 14))
# These clusters contain a variety of cells from different days:
table(cardiomyocyte_clusters@meta.data$hash.ID)
# So let's subset for only cells that are LMNA-DAY30 or WTC-DAY30
day30_cardiomyocytes <- subset(cardiomyocyte_clusters, hash.ID %in% c("LMNA-DAY30", "WTC-DAY30"))

# Get raw counts
raw_counts <- day30_cardiomyocytes@assays$RNA$counts.3


### CALCULATE CPM (TPM) =======================================================

# Calculate the TPM for each of the subsetted cells
# Calculate CPM (same as TPM for single cell as gene lengths are not needed)
cpm <- sweep(raw_counts, 2, colSums(raw_counts)/1e6, '/')

# Convert to dataframe
cpm_df <- as.data.frame(cpm)


### SAVE OUTPUT FILE ==========================================================

# save simulation output
message("Saving output to file.")
saveRDS(combined_seurat, file = snakemake@output$processed_data)
saveRDS(cpm_df, file = snakemake@output$day30_cells)


### CLEANUP ===================================================================

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)

