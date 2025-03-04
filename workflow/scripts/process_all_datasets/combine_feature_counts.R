# Script: combine_feature_counts.R

### CREATE DEBUG FILES =======================================================

# Saving image for debugging
save.image("RDA_objects/combine_feature_counts.rda")
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
  library(tidyverse)
  library(data.table)
  library(rtracklayer)
})

# Load in featureCounts results
featureCounts <- lapply(snakemake@input$featureCounts, fread, header = T, skip = 1)

# Load in metadata
metadata <- read_csv(snakemake@input$metadata)


### MATCH STAGES TO SRR ID ===================================================

# Extract Run IDs and match stages
run_ids <- sapply(snakemake@input$featureCounts, function(x) {
  str_extract(x, "SRR[0-9]+")
})

# Match stages and replicates from metadata
sample_info <- metadata %>% 
  filter(Run %in% run_ids) %>% 
  arrange(factor(Run, levels = run_ids)) %>% 
  select(Run, stage, replicate)


### CALCULATE TPM ============================================================

# Calculate TPM for each sample and organize into a matrix
tpm_matrix <- sapply(1:length(featureCounts), function(i) {
  data <- featureCounts[[i]]
  
  # Extract gene lengths and read counts
  lengths <- as.numeric(data$Length)
  reads <- as.numeric(data[[ncol(data)]])  # Last column contains read counts
  
  # Calculate reads per length
  reads_per_length <- reads / lengths
  normalization_factor <- sum(reads_per_length)
  
  # Calculate TPM for each gene
  tpm <- (reads_per_length / normalization_factor) * 1e6
  
  return(tpm)
})

# Set row and column names for the matrix
rownames(tpm_matrix) <- featureCounts[[1]]$Geneid  # Genes are in the same order across tables
colnames(tpm_matrix) <- run_ids  # Columns represent each sample


### GET THE MEAN TPM B/T REPS ================================================

# Create a list to store mean TPM values by stage
mean_tpm_by_stage <- sample_info %>%
  group_by(stage) %>%
  summarize(replicate_columns = list(Run)) %>%
  mutate(mean_tpm = map(replicate_columns, ~ {
    # Find columns in TPM matrix that correspond to the current stage's replicates
    replicate_indices <- which(colnames(tpm_matrix) %in% .x)
    # Calculate the mean TPM across replicates for each gene
    rowMeans(tpm_matrix[, replicate_indices])
  })) %>%
  select(stage, mean_tpm)

# Convert to a new data frame where each column is the mean TPM for a stage
mean_tpm_matrix <- do.call(cbind, mean_tpm_by_stage$mean_tpm)
colnames(mean_tpm_matrix) <- mean_tpm_by_stage$stage
rownames(mean_tpm_matrix) <- rownames(tpm_matrix)  # Set gene IDs as row names


### REMOVE ROWS WITH 0 TPM ACROSS ALL TIMEPOINTS ============================

# Convert mean_tpm_matrix to a data frame for filtering
mean_tpm_df <- mean_tpm_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_id")  # Convert rownames to a column for joining

# Filter out rows where all TPM values across timepoints are zero
mean_tpm_df <- mean_tpm_df %>%
  filter(rowSums(select(., -gene_id)) > 0)  # Retain only rows with non-zero TPM in any column


### ANNOTATE GENES IN MEAN TPM MATRIX ========================================

# Load annotation file (assuming it's in GTF or GFF format and contains the required fields)
annot <- import(snakemake@input$annot)

# Extract relevant columns from the annotation file and ensure unique gene_id entries
annot_df <- annot %>%
  as.data.frame() %>%
  select(gene_id, type, gene_type, gene_name) %>%
  filter(type == "gene")  # Only keep entries where type is "gene"

# Join annotation with mean TPM data
mean_tpm_annotated <- mean_tpm_df %>%
  left_join(annot_df, by = "gene_id")

# Summary: Number of entries where type == "gene"
num_types <- mean_tpm_annotated %>% filter(type == "gene") %>% nrow()
message("Number of entries with type 'gene': ", num_types)

# Summary: Number of entries where gene_type == "protein_coding"
num_gene_types <- mean_tpm_annotated %>% filter(gene_type == "protein_coding") %>% nrow()
message("Number of entries with type 'protein_coding': ", num_gene_types)


### FILTER FOR PROTEIN-CODING GENES ==========================================

# Filter mean_tpm_annotated for protein-coding genes only
mean_tpm_annotated_protein_coding <- mean_tpm_annotated %>%
  filter(gene_type == "protein_coding")

# Convert to a matrix for protein-coding genes only
mean_tpm_matrix_protein_coding <- mean_tpm_annotated_protein_coding %>%
  select(-gene_id, -type, -gene_type, -gene_name) %>%  # Exclude non-TPM columns
  as.matrix()

# Set rownames to gene IDs for mean_tpm_matrix_protein_coding
rownames(mean_tpm_matrix_protein_coding) <- mean_tpm_annotated_protein_coding$gene_id


### PLOT TPM DISTRIBUTIONS FOR PROTEIN-CODING GENES ACROSS STAGES ============

# Convert protein-coding mean TPM matrix to long format for plotting
mean_tpm_df_protein_coding <- mean_tpm_matrix_protein_coding %>%
  as.data.frame() %>%
  rownames_to_column(var = "GeneID") %>%
  pivot_longer(cols = -GeneID, names_to = "stage", values_to = "mean_tpm")

# Filter data for "day 90" stage for protein-coding genes
stage_data_protein_coding <- mean_tpm_df_protein_coding %>%
  filter(stage == "day 90")

# Ensure that mean_tpm is numeric
stage_data_protein_coding$mean_tpm <- as.numeric(stage_data_protein_coding$mean_tpm)

# Calculate thresholds for different gene numbers
gene_numbers <- c(12000, 13000, 14000, 15000, 16000)
thresholds <- sapply(gene_numbers, function(n) {
  stage_data_protein_coding %>%
    arrange(desc(mean_tpm)) %>%
    slice_head(n = n) %>%
    pull(mean_tpm) %>%
    min()
})

# Create the plot
top_n_genes_w_tpm <- ggplot(stage_data_protein_coding, aes(x = mean_tpm)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  scale_x_log10() +
  # Add vertical lines for each threshold
  geom_vline(xintercept = thresholds[1], color = "red", linetype = "dashed") +
  geom_vline(xintercept = thresholds[2], color = "orange", linetype = "dashed") +
  geom_vline(xintercept = thresholds[3], color = "darkgreen", linetype = "dashed") +
  geom_vline(xintercept = thresholds[4], color = "blue", linetype = "dashed") +
  geom_vline(xintercept = thresholds[5], color = "purple", linetype = "dashed") +
  # Add labels for each threshold
  geom_text(aes(x = thresholds[1], y = 0, label = paste("12k:", round(thresholds[1], 2))), 
            color = "red", angle = 90, vjust = -0.5, hjust = -0.1) +
  geom_text(aes(x = thresholds[2], y = 0, label = paste("13k:", round(thresholds[2], 2))), 
            color = "orange", angle = 90, vjust = -0.5, hjust = -0.1) +
  geom_text(aes(x = thresholds[3], y = 0, label = paste("14k:", round(thresholds[3], 2))), 
            color = "darkgreen", angle = 90, vjust = -0.5, hjust = -0.1) +
  geom_text(aes(x = thresholds[4], y = 0, label = paste("15k:", round(thresholds[4], 2))), 
            color = "blue", angle = 90, vjust = -0.5, hjust = -0.1) +
  geom_text(aes(x = thresholds[5], y = 0, label = paste("16k:", round(thresholds[5], 2))), 
            color = "purple", angle = 90, vjust = -0.5, hjust = -0.1) +
  labs(
    title = "TPM Distribution for Stage: Day 90 (Protein-Coding Genes n = ~19k)",
    x = "Mean TPM (log scale)",
    y = "Density"
  ) +
  theme_classic()

# Save the plot
ggsave(
  plot = top_n_genes_w_tpm,
  filename = snakemake@output$top_n_genes_w_tpm,
  device = "pdf",
  height = 5,
  width = 6.5
)


### SUBSET FOR DAYS 9, 14, 30, AND 90 =======================================

# Filter mean_tpm_annotated_protein_coding for only days 9, 14, 30, and 90
mean_tpm_subset <- mean_tpm_annotated_protein_coding %>%
  select(gene_id, gene_name, gene_type, `day 9`, `day 14`, `day 30`, `day 90`)  # Select only required days

# Remove genes with 0 TPM across all four days
mean_tpm_subset <- mean_tpm_subset %>%
  filter(rowSums(select(., `day 9`, `day 14`, `day 30`, `day 90`)) > 0)

# Count the remaining genes
remaining_genes_count <- nrow(mean_tpm_subset)
message("Number of genes with non-zero TPM in at least one of days 9, 14, 30, or 90: ", remaining_genes_count)


### PLOT TPM DISTRIBUTION ACROSS SELECTED STAGES ============================

# Convert to long format for plotting
mean_tpm_long <- mean_tpm_subset %>%
  pivot_longer(cols = starts_with("day"), names_to = "stage", values_to = "TPM") %>%
  mutate(stage = factor(stage, levels = c("day 9", "day 14", "day 30", "day 90")))

tpm_distr_per_day <- ggplot(mean_tpm_long, aes(x = TPM, fill = stage)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  labs(
    title = "TPM Distribution for Days 9, 14, 30, and 90 (Protein-Coding Genes)",
    x = "TPM (log scale)",
    y = "Density"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

# Save the plot
ggsave(
  plot = tpm_distr_per_day,
  filename = snakemake@output$tpm_distr_per_day,
  device = "pdf",
  height = 5,
  width = 6
)


### PLOT VARIANCE IN EXPRESSION FOR EACH DAY =================================

# First, get the TPM threshold from day 90 for top 14k genes
day90_threshold <- mean_tpm_subset %>%
  arrange(desc(`day 90`)) %>%
  slice_head(n = 14000) %>%
  pull(`day 90`) %>%
  min()

# Calculate genes below threshold for each day among the top 14k from day 90
top_genes_day90 <- mean_tpm_subset %>%
  arrange(desc(`day 90`)) %>%
  slice_head(n = 14000) %>%
  pull(gene_id)

genes_below_threshold <- mean_tpm_subset %>%
  filter(gene_id %in% top_genes_day90) %>%
  summarise(
    day9_below = sum(`day 9` < day90_threshold),
    day14_below = sum(`day 14` < day90_threshold),
    day30_below = sum(`day 30` < day90_threshold)
  )

# Create summary data with updated coloring logic
expression_summary <- mean_tpm_subset %>%
  arrange(desc(`day 90`)) %>%
  mutate(rank = row_number()) %>%
  # Calculate variance from day 90 values
  mutate(
    ref_expr = `day 90`,  # Use day 90 as reference
    sd_from_90 = apply(select(., `day 9`, `day 14`, `day 30`), 1, function(x) sd(c(x))),
    upper = ref_expr + sd_from_90,
    lower = ref_expr - sd_from_90
  ) %>%
  # Log transform
  mutate(
    across(c(ref_expr, sd_from_90, upper, lower, `day 9`, `day 14`, `day 30`, `day 90`), 
           ~log10(.x + 1))
  ) %>%
  # Convert to long format for plotting individual days
  pivot_longer(
    cols = starts_with("day"),
    names_to = "stage",
    values_to = "TPM"
  ) %>%
  # Add column to determine if point should be colored
  mutate(
    point_color = case_when(
      stage == "day 90" ~ "gray80",  # day 90 always gray
      rank <= 14000 & TPM < log10(day90_threshold + 1) ~ stage,  # top ranked genes: color if below threshold
      rank > 14000 & TPM > log10(day90_threshold + 1) ~ stage,   # lower ranked genes: color if above threshold
      TRUE ~ "gray80"  # everything else gray
    )
  )

# Create the plot
expression_variance <- ggplot(expression_summary, aes(x = rank)) +
  # Add shaded area for variance from day 90
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.01, fill = "gray50") +
  # Add points with conditional coloring
  geom_point(aes(y = TPM, color = point_color), size = 0.1) +
  # Add day 90 reference line
  geom_line(aes(y = ref_expr), color = "black", size = 0.8) +
  # Add vertical line at 14k genes
  geom_vline(xintercept = 14000, linetype = "dashed", color = "red", alpha = 0.5) +
  # Add horizontal line for day 90 threshold
  geom_hline(yintercept = log10(day90_threshold + 1), linetype = "dashed", color = "blue", alpha = 0.5) +
  # Add annotation for threshold
  annotate("text", x = max(expression_summary$rank), 
           y = log10(day90_threshold + 1), 
           label = paste("Day 90 threshold TPM:", round(day90_threshold, 2)), 
           hjust = 1.8, vjust = -5.5, color = "blue") +
  scale_color_manual(
    values = c(
      "gray80" = "gray80",
      "day 9" = "#E41A1C",
      "day 14" = "#377EB8",
      "day 30" = "#4DAF4A"
    )
  ) +
  labs(
    title = "Gene Expression Profile Across Days",
    subtitle = paste0("Ranked by Day 90 expression (black line)\n",
                      "For top 14k genes: colored if below threshold\n",
                      "For remaining genes: colored if above threshold"),
    x = "Gene Rank",
    y = "log10(TPM + 1)",
    color = "Stage"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

# Save the plot
ggsave(
  plot = expression_variance,
  filename = snakemake@output$expression_variance,
  device = "pdf",
  height = 5,
  width = 5
)


### VISUALIZE GENE DIFFERENCES ================================================

# Count the number of pairs that are removed on each day according to this cutoff
threshold_comparison <- mean_tpm_subset %>%
  filter(gene_id %in% top_genes_day90) %>%
  summarise(
    day9_below = sum(`day 9` < day90_threshold),
    day14_below = sum(`day 14` < day90_threshold),
    day30_below = sum(`day 30` < day90_threshold)
  )

# Create data frame for plotting
plot_data <- data.frame(
  Day = factor(c("Day 9", "Day 14", "Day 30"), levels = c("Day 9", "Day 14", "Day 30")),
  Count = c(threshold_comparison$day9_below, 
            threshold_comparison$day14_below, 
            threshold_comparison$day30_below)
)

# Create the plot
threshold_plot <- ggplot(plot_data, aes(x = Day, y = Count, fill = Day)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5) +
  scale_fill_manual(values = c("#E41A1C", "#377EB8", "#4DAF4A")) +
  labs(
    title = "Genes Below Day 90 Threshold",
    subtitle = paste("Number of day 90 top 14k genes with TPM <", round(day90_threshold, 2), "in each day"),
    x = "Time Point",
    y = "Number of Genes"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11)
  )

# Save the plot
ggsave(
  plot = threshold_plot,
  filename = snakemake@output$threshold_plot,
  device = "pdf",
  height = 4,
  width = 5
)


### CONFIRMING CARDIOMYOCYTE GENE PRESENCE IN GENE LIST ======================

# Load in the cardiomyocyte gene list from PangaoDB
cardiomyocyte_genes <- read_tsv(snakemake@input$cardiomyocyte_genes) %>% 
  filter(`cell type` == "Cardiomyocytes")

# Get gene symbols
cardio_gene_list <- cardiomyocyte_genes$`official gene symbol`
# Check presence of cardio_gene_list in mean_tpm_subset$gene_name
message(paste0("Genes in cardiomyocyte gene list that don't have any reads on days 9, 14, 30, or 90: ", paste(cardio_gene_list[!cardio_gene_list %in% mean_tpm_subset$gene_name], collapse=", ")))
# Genes in cardiomyocyte gene list that don't have any reads on days 9, 14, 30, or 90: COX8B, MHRT
cardio_gene_list <- cardio_gene_list[cardio_gene_list %in% mean_tpm_subset$gene_name]

# Do all of these genes pass the TPM filter when subsetting for top 14k genes in day 90?
tpm_genes_in_cardio_gene_list <- mean_tpm_subset %>%
  filter(gene_name %in% cardio_gene_list)
  
# Create a long-format table with a "pass_threshold" column to indicate if TPM > threshold for each gene/day
tpm_genes_long <- tpm_genes_in_cardio_gene_list %>%
  pivot_longer(cols = c(`day 9`, `day 14`, `day 30`, `day 90`), names_to = "day", values_to = "TPM") %>%
  mutate(log_TPM = log10(TPM + 1),  # Log-transform TPM to avoid log(0) issues
         pass_threshold = TPM > day90_threshold)

# Filter to include only genes that fail the threshold in at least one day
genes_failing_threshold <- tpm_genes_long %>%
  group_by(gene_name) %>%
  filter(any(!pass_threshold)) %>%  # Keep only genes with at least one "false" in pass_threshold
  ungroup()

# Plot only genes that fail the threshold at least once
failed_thresh_at_least_once <- ggplot(genes_failing_threshold, aes(x = gene_name, y = TPM, color = pass_threshold)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = day90_threshold, linetype = "dashed", color = "red") +
  facet_wrap(~ day, scales = "free_x") +  # Separate plot for each day
  labs(
    title = "TPM for Cardiomyocyte Genes Failing Threshold in at Least One Day",
    x = "Gene Name",
    y = "TPM"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save the plot
ggsave(
  plot = failed_thresh_at_least_once,
  filename = snakemake@output$failed_thresh_at_least_once,
  device = "pdf",
  height = 5,
  width = 5
)

genes_in_all_days <- setdiff(tpm_genes_in_cardio_gene_list$gene_name, unique(genes_failing_threshold$gene_name))
genes_in_no_days <- unique((tpm_genes_long %>% group_by(gene_name) %>% filter(all(!pass_threshold)) %>% ungroup())$gene_name)
genes_not_in_day_90 <- unique(tpm_genes_long %>% filter(day == "day 90", pass_threshold == FALSE) %>% pull(gene_name))

# Print useful information:
message(paste0("Genes in the cardiomyocyte gene set that are present all days: ", paste(genes_in_all_days, collapse = ", ")))
message(paste0("Genes in the cardiomyocyte gene set that are absent all days: ", paste(genes_in_no_days, collapse = ", ")))
message(paste0("Genes in the cardiomyocyte gene set that are absent on day 90: ", paste(genes_not_in_day_90, collapse = ", ")))


### CARDIOMYOCYTE GENES THAT DIFFER BETWEEN DAYS ==============================

# Filter for only cardiomyocyte genes
cardio_tpm_data <- mean_tpm_subset %>%
  filter(gene_name %in% cardio_gene_list)

# Find genes that:
# 1. Don't meet threshold on day 90
# 2. But meet threshold in at least one earlier day
cardio_genes_analysis <- cardio_tpm_data %>%
  filter(`day 90` < day90_threshold) %>%  # Don't meet threshold on day 90
  filter(`day 9` >= day90_threshold | 
           `day 14` >= day90_threshold | 
           `day 30` >= day90_threshold)  # Meet threshold in at least one early day

# Create summary
gene_summary <- cardio_genes_analysis %>%
  select(gene_name, `day 9`, `day 14`, `day 30`, `day 90`)

# Print results
if (nrow(gene_summary) > 0) {
  print("Cardiomyocyte genes that meet day 90 threshold in earlier days but not on day 90:")
  print(gene_summary)
} else {
  print("No cardiomyocyte genes meet the criteria")
}

# Create a longer format for plotting
gene_summary_long <- gene_summary %>%
  pivot_longer(
    cols = starts_with("day"),
    names_to = "day",
    values_to = "TPM"
  ) %>%
  mutate(day = factor(day, levels = c("day 9", "day 14", "day 30", "day 90")))  # Order the days

# Create visualization
cardio_genes_that_decrease_expr <- ggplot(gene_summary_long, aes(x = day, y = TPM, group = gene_name, color = gene_name)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_hline(yintercept = day90_threshold, linetype = "dashed", color = "red") +
  labs(
    title = "Cardiomyocyte Genes Meeting day 90 Threshold Earlier but not on Day 90",
    subtitle = paste("Day 90 threshold:", round(day90_threshold, 2)),
    x = "Day",
    y = "TPM",
    color = "Gene"
  ) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save the plot
ggsave(
  plot = cardio_genes_that_decrease_expr,
  filename = snakemake@output$cardio_genes_that_decrease_expr,
  device = "pdf",
  height = 4,
  width = 5
)


### CREATING OUTPUT ===========================================================

# Define the days and thresholds
days <- c("day 9", "day 14", "day 30", "day 90")
thresholds <- c(12000, 13000, 14000, 15000, 16000)

# Initialize a list to store the thresholds for each day
threshold_values <- list()

# Calculate TPM thresholds for each day and each threshold value
for (day in days) {
  tpm_values <- mean_tpm_subset[[day]]
  tpm_values_sorted <- sort(tpm_values, decreasing = TRUE)
  
  for (thresh in thresholds) {
    if (length(tpm_values_sorted) >= thresh) {
      threshold_tpm <- tpm_values_sorted[thresh]
    } else {
      threshold_tpm <- min(tpm_values_sorted)
    }
    # Store the threshold TPM value
    threshold_values[[paste0(day, "_top", thresh)]] <- threshold_tpm
  }
}

# Create TRUE/FALSE columns for each day and threshold
for (day in days) {
  tpm_values <- mean_tpm_subset[[day]]
  
  for (thresh in thresholds) {
    threshold_tpm <- threshold_values[[paste0(day, "_top", thresh)]]
    # Create a new column name, e.g., 'day 9_top13000'
    col_name <- paste0(day, "_top", thresh)
    mean_tpm_subset[[col_name]] <- tpm_values >= threshold_tpm
  }
}

# Select the desired columns for the final output dataframe
columns_to_select <- c("gene_id", "gene_name", days,
                       paste0(rep(days, each = length(thresholds)), "_top", thresholds))

# Create the final dataframe
df <- mean_tpm_subset[, columns_to_select]


## SAVE OUTPUT FILE =========================================================

# save simulation output
message("Saving output to file.")
saveRDS(df, file = snakemake@output$gene_list_df)


### CLEANUP =================================================================

# close log file connection
message("Closing log file connection")
sink()
sink(type = "message")
close(log)
