# Script: processing_raw_counts_data.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/processing_raw_counts_data.rda"))
message("Saved Image")
# stop("Manually Stopped Program after Saving Image")

# Open log file to collect messages, warnings, and errors
log_filename <- snakemake@log[[1]]
log <- file(log_filename, open = "wt")
sink(log)
sink(log, type = "message")


### LOADING FILES =============================================================

message("Loading in packages")
suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(SoupX)
  library(SeuratWrappers)
})

message("Loading input files")
umi_files <- snakemake@input$umi_files
raw_matrices <- snakemake@input$raw_matrices
filtered_matrices <- snakemake@input$filt_matrices

### QCs ======================================================================

for (f in umi_files) {
  
  # Extract sample name, e.g. "P1_24s001312-1-1_Xie_lane124s001312"
  sample_name <- basename(dirname(dirname(dirname(f))))
  
  # Read UMIperCell.txt (2 cols: UMI_count, is_filtered)
  df <- read_tsv(f, col_names = c("UMI_count", "is_filtered"))
  
  # Because the file is already sorted by descending UMIs, define the rank:
  df$rank <- seq_len(nrow(df))
  
  # Plot with log-scale on both x and y, coloring by is_filtered
  p <- ggplot(df, aes(x = rank, y = UMI_count, color = factor(is_filtered))) +
    geom_line() +
    scale_x_log10() + 
    scale_y_log10() +
    scale_color_manual(
      values = c("0" = "gray60", "1" = "blue"),
      labels = c("Background", "Filtered Cell")
    ) +
    labs(
      title = "Barcode Rank (Knee) Plot",
      subtitle = sample_name,
      x = "Barcode Rank (log10)",
      y = "UMI Count (log10)",
      color = NULL
    ) +
    theme_classic()
  
  # Save or print the plot
  ggsave(filename = file.path(dirname(f), "barcode_rank_plot.pdf"),
         plot = p, device = "pdf", width = 5, height = 4)
}


### READ IN DATA =============================================================

raw_data_list <- lapply(dirname(raw_matrices), Read10X)
names(raw_data_list) <- basename(dirname(dirname(dirname(dirname(raw_matrices)))))

filtered_data_list <- lapply(dirname(filtered_matrices), Read10X)
names(filtered_data_list) <- basename(dirname(dirname(dirname(dirname(filtered_matrices)))))


### CREATE SEURAT OBJECTS ====================================================

# Create SoupChannel objects
soup_channel_list <- lapply(names(raw_data_list), function(name) {
  filtered_name <- gsub("raw", "filtered", name)
  raw_data <- raw_data_list[[name]]
  filtered_data <- filtered_data_list[[filtered_name]]
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
  output_dir <- dirname(snakemake@output$combined_seurat)
  
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
  
  # Create Seurat object with adjusted counts and CMO data
  seurat_obj <- CreateSeuratObject(counts = adjusted_counts)
  return(seurat_obj)
})
names(seurat_objects) <- names(adjusted_counts_list)


### QUALITY CONTROL ==========================================================

# Calculate the percentage of mitochondrial & ribosomal genes
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  # Percentage of mitochondrial genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  
  # Percentage of ribosomal genes
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
  
  return(seurat_obj)
})


### VISUALIZE QCs ============================================================

# Label each object with a sample identifier
seurat_objects[[1]]$sample_id <- "Sample1"
seurat_objects[[2]]$sample_id <- "Sample2"

# Merge the two Seurat objects into one and normalize
merged_obj <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[[2]],
  project = "MergedSamples"
)
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)

# Function to plot QC metrics as violin plots
plot_violin_qc <- function(seurat_obj, qc_metric) {
  # Ensure data slot is normalized
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  
  # Check if the QC metric exists in the metadata
  if (!(qc_metric %in% colnames(seurat_obj@meta.data))) {
    stop(paste("QC metric", qc_metric, "not found in Seurat object's metadata."))
  }
  
  # Generate violin plot
  p <- VlnPlot(
    seurat_obj,
    features = qc_metric,
    group.by = "sample_id",  # Assumes sample_id column exists for grouping
    pt.size = 0.1
  ) +
    NoLegend()
}

# Create individual QC plots before filtering
plot1_before <- plot_violin_qc(merged_obj, "nFeature_RNA")
plot2_before <- plot_violin_qc(merged_obj, "nCount_RNA")
plot3_before <- plot_violin_qc(merged_obj, "percent.mt")
plot4_before <- plot_violin_qc(merged_obj, "percent.ribo")

# Define output directory based on combined_seurat output path
output_dir <- dirname(snakemake@output$combined_seurat)

# Save each "before QC" plot as PDF
ggsave(filename = file.path(output_dir, "nFeature_RNA_before_qc.pdf"),
       plot = plot1_before, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "nCount_RNA_before_qc.pdf"),
       plot = plot2_before, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "percent.mt_before_qc.pdf"),
       plot = plot3_before, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "percent.ribo_before_qc.pdf"),
       plot = plot4_before, device = "pdf", width = 5, height = 4)


### QUALITY CONTROL ==========================================================

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

# Label each object with a sample identifier
seurat_objects[[1]]$sample_id <- "Sample1"
seurat_objects[[2]]$sample_id <- "Sample2"

# Merge the two Seurat objects into one and normalize
merged_obj <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[[2]],
  project = "MergedSamples"
)
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)

# Modified to include "(After QC)"
plot_violin_qc_after <- function(seurat_obj, qc_metric) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  if (!(qc_metric %in% colnames(seurat_obj@meta.data))) {
    stop(paste("QC metric", qc_metric, "not found in Seurat object's metadata."))
  }
  
  p <- VlnPlot(
    seurat_obj,
    features = qc_metric,
    group.by = "sample_id",
    pt.size = 0.1
  ) + NoLegend() + theme_classic() +
    ggtitle(paste(qc_metric, "(After QC)"))
  
  return(p)
}

plot1_after <- plot_violin_qc_after(merged_obj, "nFeature_RNA")
plot2_after <- plot_violin_qc_after(merged_obj, "nCount_RNA")
plot3_after <- plot_violin_qc_after(merged_obj, "percent.mt")
plot4_after <- plot_violin_qc_after(merged_obj, "percent.ribo")

# Save each "after QC" plot as PDF
ggsave(filename = file.path(output_dir, "nFeature_RNA_after_qc.pdf"),
       plot = plot1_after, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "nCount_RNA_after_qc.pdf"),
       plot = plot2_after, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "percent.mt_after_qc.pdf"),
       plot = plot3_after, device = "pdf", width = 5, height = 4)
ggsave(filename = file.path(output_dir, "percent.ribo_after_qc.pdf"),
       plot = plot4_after, device = "pdf", width = 5, height = 4)


### FILTER MITOCHONDRIAL CELLS ================================================

# Filter out cells with >50% mitochondrial reads
seurat_objects <- lapply(seurat_objects, function(seurat_obj) {
  subset(seurat_obj, subset = percent.mt <= 50)
})

# How many cells are left
message("Number of cells left after filtering out likely dead cells")
lapply(seurat_objects, ncol)


### MERGE AND NORMALIZE RUNS =================================================

# Merge the Seurat objects into one combined object
print(seurat_objects)
combined_seurat <- merge(
  x = seurat_objects[[1]],
  y = seurat_objects[[2]],
  project = "Combined_Seurat"
)
print(combined_seurat)

# Normalize the RNA assay counts
combined_seurat <- NormalizeData(combined_seurat, assay = "RNA")
combined_seurat <- JoinLayers(combined_seurat, overwrite = TRUE)
combined_seurat <- ScaleData(combined_seurat)
print(combined_seurat)


### IDENTIFY VARIABLE FEATURES AND INTEGRATE BATCHES =========================

# Split the Seurat object into a list by 'run' (dataset) and identify variable features in each run
message("Identifying variable features in each sample")
seurat_list <- SplitObject(combined_seurat, split.by = "sample_id")
seurat_list <- lapply(seurat_list, FindVariableFeatures)

# Select integration features across datasets
message("Selecting integration features")
features <- SelectIntegrationFeatures(object.list = seurat_list)

# Perform integration using RunFastMNN, specifying 'individual' as the batch variable
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
umap_sample <- DimPlot(combined_seurat, reduction = "umap", group.by = "sample_id", label = TRUE) +  NoLegend() + ggtitle("UMAP: Sample IDs")

# Save UMAP plots
ggsave(filename = file.path(output_dir, "umap_clusters.pdf"),
       plot = umap_clusters, device = "pdf", width = 6, height = 5)
ggsave(filename = file.path(output_dir, "umap_sample.pdf"),
       plot = umap_sample, device = "pdf", width = 6, height = 5)


### IDENTIFYING CARDIOMYOCYTE SPECIFIC CLUSTER ===============================

# Specific markers that Liangfu sent
cardiac_tfs <- c("TNNT2", "MYH6", "MYH7", "ACTC1", "TPM1")
cardiac_feature_plot <- FeaturePlot(combined_seurat, 
                                    features = cardiac_tfs, 
                                    ncol = 3, 
                                    order = TRUE, 
                                    min.cutoff = 'q10', 
                                    max.cutoff = 'q90')

# Save the cardiac features plot
ggsave(filename = file.path(output_dir, "cardiac_feature_plot.pdf"),
       plot = cardiac_feature_plot, device = "pdf", width = 8, height = 6)


### CALCULATE CPM (TPM) =======================================================

# Extract counts.1 and counts.2
counts_1 <- combined_seurat@assays$RNA$counts.1
counts_2 <- combined_seurat@assays$RNA$counts.2
# Combine the two sparse matrices by columns (since they share the same features/rows)
combined_counts <- cbind(counts_1, counts_2)

# Calculate the TPM for each of the subsetted cells
# Calculate CPM (same as TPM for single cell as gene lengths are not needed)
cpm <- sweep(combined_counts, 2, colSums(combined_counts)/1e6, '/')

# Convert to dataframe
cpm_df <- as.data.frame(cpm)


### SAVE OUTPUT ===============================================================

message("Saving output files")
saveRDS(combined_seurat, file = snakemake@output$combined_seurat)
saveRDS(cpm_df, file = snakemake@output$cpm_df)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)