# Load necessary libraries
library(opi)
library(tidyverse)
library(glue)

# Initialize the datastore connection
datastore <- opi()

# Function to fetch epilepsy studies
fetch_disease_studies <- function(organism) {
  sql <- glue("
    SELECT *
    FROM {organism}_disease_b38_gc33.comparisons
    JOIN {organism}_disease_b38_gc33.comparison_data USING (comparison_index)
    JOIN {organism}_disease_b38_gc33.gene_annotation USING (gene_index)
    WHERE LOWER(case_disease_state) LIKE '%alzheimer%' 
    AND comparison_category = 'Disease vs. Normal'
  ")
  disease_studies <- datastore$query(sql)
  return(disease_studies)
}

# Fetch epilepsy studies for human
human_disease_studies <- fetch_disease_studies('human')

# View the first few rows of the fetched data
print(head(human_disease_studies))

# Export the fetched data to a CSV file for further analysis
write_csv(human_disease_studies, "human_alzheimer_studies_differential_expression.csv")

