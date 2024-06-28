# Load necessary libraries
library(opi)
library(tidyverse)
library(broom)

# Initialize the datastore connection
datastore <- opi()

# Function to obtain TPM and MS values and calculate correlation for each gene
calculate_correlation <- function(gene_id) {
  # Query to obtain TPM values for the gene
  tpm_sql <- stringr::str_glue("
    SELECT spl.sample_id, spl.tissue_gtex AS tissue, fpkm * tpm_scaling_factor AS tpm
    FROM gtex_b38_gc33.samples spl
    JOIN gtex_b38_gc33.gene_fpkm fpkm ON fpkm.sample_index = spl.sample_index
    WHERE fpkm.gene_index = {gene_id}
  ")
  tpm_data <- datastore$query(tpm_sql) %>% as_tibble()
  
  # Query to obtain MS values for the gene
  ms_sql <- stringr::str_glue("
    SELECT spl.sample_id, spl.tissue_gtex AS tissue, ms.value AS ms
    FROM gtex_b38_gc33.samples spl
    JOIN gtex_b38_gc33.protein_ms ms ON ms.sample_index = spl.sample_index
    WHERE ms.gene_index = {gene_id}
  ")
  ms_data <- datastore$query(ms_sql) %>% as_tibble()
  
  # Join TPM and MS data by sample_id and tissue
  combined_data <- inner_join(tpm_data, ms_data, by = c("sample_id", "tissue"))
  
  # Calculate correlation for each tissue
  correlation_results <- combined_data %>%
    group_by(tissue) %>%
    summarise(
      correlation = cor(tpm, ms, use = "complete.obs"),
      p_value = cor.test(tpm, ms, use = "complete.obs")$p.value,
      .groups = 'drop'
    ) %>%
    mutate(gene_id = gene_id)
  
  return(correlation_results)
}

# List of genes to analyze
gene_list_sql <- "SELECT DISTINCT gene_index FROM gtex_b38_gc33.gene_fpkm"
gene_list <- datastore$query(gene_list_sql) %>% pull(gene_index)

# Apply the function to each gene and store results in a list
results_list <- lapply(gene_list, calculate_correlation)

# Combine all results into a single dataframe
all_results <- bind_rows(results_list)

# View the combined results
print(all_results)

# Optionally, write the combined results to a CSV file for further analysis or sharing
write_csv(all_results, "mRNA_protein_correlation_results.csv")
