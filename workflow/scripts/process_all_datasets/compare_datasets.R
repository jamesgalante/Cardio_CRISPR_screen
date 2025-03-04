# Script: # Script: compare_datasets.R

### SETUP =====================================================================

# Saving image for debugging
save.image(paste0("RDA_objects/compare_datasets.rda"))
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
  library(readxl)
})

message("Loading input files")

# Read in each dataset
sean_m_wu <- readRDS(snakemake@input$sean_m_wu)
joseph_c_wu <- readRDS(snakemake@input$joseph_c_wu)
liangfu <- readRDS(snakemake@input$liangfu)

# Read in the cardiomyocyte genes
cardiomyocyte_genes <- read_tsv(snakemake@input$cardiomyocyte_genes) %>% 
  filter(`cell type` == "Cardiomyocytes") %>% 
  pull(`official gene symbol`)

# Read in the kinase genes
kinase_genes <- read_xlsx(snakemake@input$kinase_genes) %>% pull(`Target Gene Symbol`) %>% unique()
kinase_genes <- kinase_genes[kinase_genes != "Non-Targeting Control"]

# Read in the membrane proteins
membrane_genes <- read_xlsx(snakemake@input$membrane_genes, sheet = 2, col_names = F)
membrane_genes <- membrane_genes %>%
  # Rename the columns for convenience
  dplyr::rename(label = 1, sequence = 2) %>%
  # Remove negative controls (those starting with "0")
  filter(!startsWith(label, "0")) %>%
  # Split the label into components: GeneID, Symbol, Sublibrary, GuideID
  separate(label, 
           into = c("GeneID", "Symbol", "Sublibrary", "GuideID"), 
           sep = "_", remove = FALSE) %>%
  # Keep only rows whose Sublibrary is "MEPR"
  filter(Sublibrary == "MEPR") %>%
  # Get unique symbols
  distinct(Symbol) %>%
  pull(Symbol)


### FORMATTING ================================================================

# Calculate average CPM for sean_m_wu across all cells
sean_m_wu_mean <- rowMeans(sean_m_wu)

# Calculate average CPM for liangfu dataset across all cells
liangfu_mean <- rowMeans(liangfu)

# Because sean_m_wu is day 30, let's compare with day 30 from the joseph c wu dataset
joseph_c_wu_day30 <- joseph_c_wu %>% select(gene_name, `day 30`)

# (Optional) message about the length of each:
message(paste("Number of genes in Sean M. Wu:", length(sean_m_wu_mean)))
message(paste("Number of genes in Liangfu:", length(liangfu_mean)))
message(paste("Number of genes in Joseph C. Wu Day 30:", nrow(joseph_c_wu_day30)))


### SUBSET TO COMPARE SIMILAR GENES ===========================================

# 1) Find intersect among Sean, Joseph (day30), Liangfu
common_genes_3way <- Reduce(intersect, list(
  names(sean_m_wu_mean),
  joseph_c_wu_day30$gene_name,
  names(liangfu_mean)
))
message(paste0("There are ", length(common_genes_3way), " genes in common among Sean, Joseph (day 30), and Liangfu."))

# 2) Subset each dataset to common genes
sean_subset <- sean_m_wu_mean[common_genes_3way]
liangfu_subset <- liangfu_mean[common_genes_3way]
joe_subset <- joseph_c_wu_day30 %>% 
  filter(gene_name %in% common_genes_3way)

# 3) Combine them into one table
sean_subset_df <- data.frame(gene_name = names(sean_subset), 
                             sean_day_30 = sean_subset)
liangfu_subset_df <- data.frame(gene_name = names(liangfu_subset), 
                                liangfu_day_30 = liangfu_subset)

combined_data_all_three <- joe_subset %>%
  dplyr::rename(joe_day_30 = `day 30`) %>%
  left_join(sean_subset_df, by = "gene_name") %>%
  left_join(liangfu_subset_df, by = "gene_name")


### COMPARE TPMs ==============================================================

# 1) Joe vs. Sean
plot_joe_sean <- combined_data_all_three %>%
  ggplot(aes(x = joe_day_30, y = sean_day_30)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  labs(title = "Joe vs. Sean (Day 30)", x = "Joseph C. Wu Day 30", y = "Sean M. Wu Day 30")
ggsave(plot = plot_joe_sean, filename = "results/compare_datasets/plot_joe_sean.pdf", device = "pdf", height = 4, width = 4)

# 2) Joe vs. Liangfu
plot_joe_liangfu <- combined_data_all_three %>%
  ggplot(aes(x = joe_day_30, y = liangfu_day_30)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  labs(title = "Joe vs. Liangfu (Day 30)", x = "Joseph C. Wu Day 30", y = "Liangfu Day 30")
ggsave(plot = plot_joe_liangfu, filename = "results/compare_datasets/plot_joe_liangfu.pdf", device = "pdf", height = 4, width = 4)


# 3) Sean vs. Liangfu
plot_sean_liangfu <- combined_data_all_three %>%
  ggplot(aes(x = sean_day_30, y = liangfu_day_30)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = 'lm', formula = y ~ x) +
  scale_x_log10() + scale_y_log10() +
  theme_classic() +
  labs(title = "Sean vs. Liangfu (Day 30)", x = "Sean M. Wu Day 30", y = "Liangfu Day 30")
ggsave(plot = plot_sean_liangfu, filename = "results/compare_datasets/plot_sean_liangfu.pdf", device = "pdf", height = 4, width = 4)



### CALCULATE CUTOFF SIMILARITY ===============================================

calculate_3way_overlap <- function(data, n_genes) {
  data <- data %>%
    mutate(
      in_top_joe = rank(-joe_day_30) <= n_genes,
      in_top_sean = rank(-sean_day_30) <= n_genes,
      in_top_liangfu = rank(-liangfu_day_30) <= n_genes
    )
  
  # 3-way intersection
  all_three <- sum(data$in_top_joe & data$in_top_sean & data$in_top_liangfu)
  
  # Only Joe
  only_joe <- sum(data$in_top_joe & !data$in_top_sean & !data$in_top_liangfu)
  
  # Only Sean
  only_sean <- sum(!data$in_top_joe & data$in_top_sean & !data$in_top_liangfu)
  
  # Only Liangfu
  only_liangfu <- sum(!data$in_top_joe & !data$in_top_sean & data$in_top_liangfu)
  
  # Joe & Sean but not Liangfu
  joe_and_sean <- sum(data$in_top_joe & data$in_top_sean & !data$in_top_liangfu)
  
  # Joe & Liangfu but not Sean
  joe_and_liangfu <- sum(data$in_top_joe & !data$in_top_sean & data$in_top_liangfu)
  
  # Sean & Liangfu but not Joe
  sean_and_liangfu <- sum(!data$in_top_joe & data$in_top_sean & data$in_top_liangfu)
  
  # Return a named vector or a small data frame
  overlap_stats <- c(
    n_genes = n_genes,
    all_three = all_three,
    only_joe = only_joe,
    only_sean = only_sean,
    only_liangfu = only_liangfu,
    joe_and_sean = joe_and_sean,
    joe_and_liangfu = joe_and_liangfu,
    sean_and_liangfu = sean_and_liangfu
  )
  return(overlap_stats)
}

thresholds <- c(13000, 14000, 15000, 16000)
for (n in thresholds) {
  stats_3way <- calculate_3way_overlap(combined_data_all_three, n)
  message(sprintf("\nFor top %d genes (3-way comparison):", n))
  print(stats_3way)
}


### CARDIOMYOCYTE SPECIFIC GENES ? ============================================

calculate_overlap_with_cardio_3way <- function(data, n_genes, cardio_genes) {
  data <- data %>%
    mutate(
      in_top_joe = rank(-joe_day_30) <= n_genes,
      in_top_sean = rank(-sean_day_30) <= n_genes,
      in_top_liangfu = rank(-liangfu_day_30) <= n_genes,
      is_cardio = gene_name %in% cardio_genes
    )
  
  # Only in Joe's top N
  only_joe_cardio <- data %>%
    filter(in_top_joe & !in_top_sean & !in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  # Only in Sean's top N
  only_sean_cardio <- data %>%
    filter(!in_top_joe & in_top_sean & !in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  # Only in Liangfu's top N
  only_liangfu_cardio <- data %>%
    filter(!in_top_joe & !in_top_sean & in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  # In all three
  all_three_cardio <- data %>%
    filter(in_top_joe & in_top_sean & in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  # (Optional) In two out of three, e.g. Joe & Sean only, etc.
  # E.g. Joe & Sean but not Liangfu:
  joe_and_sean_cardio <- data %>%
    filter(in_top_joe & in_top_sean & !in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  # None of the three
  none_cardio <- data %>%
    filter(!in_top_joe & !in_top_sean & !in_top_liangfu & is_cardio) %>%
    pull(gene_name)
  
  message(sprintf("\nFor top %d genes (3-way with cardiomyocyte genes):", n_genes))
  
  message(paste("Cardio genes only in Joe's top N:", ifelse(length(only_joe_cardio) == 0, "None", paste(only_joe_cardio, collapse=", "))))
  message(paste("Cardio genes only in Sean's top N:", ifelse(length(only_sean_cardio) == 0, "None", paste(only_sean_cardio, collapse=", "))))
  message(paste("Cardio genes only in Liangfu's top N:", ifelse(length(only_liangfu_cardio) == 0, "None", paste(only_liangfu_cardio, collapse=", "))))
  # message(paste("Cardio genes in top N of all three:", ifelse(length(all_three_cardio) == 0, "None", paste(all_three_cardio, collapse=", "))))
  message(paste("Cardio genes in none of the three datasets' top N:", ifelse(length(none_cardio) == 0, "None", paste(none_cardio, collapse=", "))))
}

# Run function to take the top N genes
thresholds <- c(13000, 14000, 15000, 16000)
for(n in thresholds) {
  calculate_overlap_with_cardio_3way(combined_data_all_three, n, cardiomyocyte_genes)
}


### CREATE FINAL GENE LIST ====================================================

# Subset for the top 13k
final_top_13k <- combined_data_all_three %>%
  mutate(
    in_top_joe = rank(-joe_day_30) <= 13000,
    in_top_sean = rank(-sean_day_30) <= 13000,
    in_top_liangfu = rank(-liangfu_day_30) <= 13000,
    is_cardio = gene_name %in% cardiomyocyte_genes
  ) %>%
  filter(in_top_joe, in_top_sean, in_top_liangfu) %>%
  pull(gene_name)


### ADD EXTRA DATASETS ========================================================

# Union the cardiomyocyte genes and the final gene list to include all those genes even if they didn't make the TPM cutoff
final_top_13k <- union(cardiomyocyte_genes, final_top_13k)

# Union the kinase genes and the final gene list
final_top_13k <- union(kinase_genes, final_top_13k)

# Union the membrane genes and the final gene list
final_top_13k <- union(membrane_genes, final_top_13k)


### SAVE OUTPUT ===============================================================

# Save output files
message("Saving output files")
writeLines(final_top_13k, snakemake@output$final_top_13k)


### CLEAN UP ==================================================================

message("Closing log file")
sink()
sink(type = "message")
close(log)