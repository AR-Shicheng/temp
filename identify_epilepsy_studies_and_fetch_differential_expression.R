# Load necessary libraries
library(opi)
library(tidyverse)
library(glue)

# Initialize the datastore connection
datastore <- opi()

# Function to fetch epilepsy studies
fetch_epilepsy_studies <- function(organism) {
  sql <- glue("
    SELECT *
    FROM {organism}_disease_b6.comparisons
    JOIN {organism}_disease_b6.comparison_data USING (comparison_index)
    JOIN {organism}_disease_b6.gene_annotation USING (gene_index)
    WHERE case_disease_state LIKE '%epilep%'
    LIMIT 10
  ")
  epilepsy_studies <- datastore$query(sql)
  return(epilepsy_studies)
}

# Fetch epilepsy studies for rat
rat_epilepsy_studies <- fetch_epilepsy_studies('rat')

# View the first few rows of the fetched data
print(head(rat_epilepsy_studies))

# Export the fetched data to a CSV file for further analysis
write_csv(rat_epilepsy_studies, "rat_epilepsy_studies_differential_expression.csv")

# Function to fetch and join differential expression data
fetch_differential_expression <- function(organism, disease_state) {
  sql <- glue("
    SELECT
      comparison_data.gene_index,
      gene_annotation.gene_name,
      comparison_data.log2fc,
      comparison_data.p_value,
      comparisons.case_disease_state,
      comparisons.control_disease_state,
      comparisons.comparison_index
    FROM {organism}_disease_b6.comparisons AS comparisons
    JOIN {organism}_disease_b6.comparison_data AS comparison_data
    ON comparisons.comparison_index = comparison_data.comparison_index
    JOIN {organism}_disease_b6.gene_annotation AS gene_annotation
    ON comparison_data.gene_index = gene_annotation.gene_index
    WHERE comparisons.case_disease_state LIKE '%{disease_state}%'
  ")
  
  differential_expression_data <- datastore$query(sql)
  return(differential_expression_data)
}

# Fetch differential expression data for epilepsy in rats
rat_epilepsy_diff_expression <- fetch_differential_expression('rat', 'epilep')

# View the first few rows of the differential expression data
print(head(rat_epilepsy_diff_expression))

# Export the differential expression data to a CSV file for further analysis
write_csv(rat_epilepsy_diff_expression, "rat_epilepsy_differential_expression.csv")
