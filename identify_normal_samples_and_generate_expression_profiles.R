# Load necessary libraries
library(opi)
library(tidyverse)

# Initialize the datastore connection
datastore <- opi()

# Function to fetch control samples
fetch_control_samples <- function(organism) {
  sql <- glue::glue("
    SELECT *
    FROM {organism}_disease_b6.comparisons
    WHERE comparisons.organism = '{organism}'
      AND comparisons.comparison_category = 'Disease vs. Normal'
      AND comparisons.control_subject_treatment IS NULL
      AND comparisons.control_cell_type IS NULL
    LIMIT 10
  ")
  control_samples <- datastore$query(sql)
  return(control_samples)
}

# Fetch control samples for rat, mouse, and human
rat_control_samples <- fetch_control_samples('rat')
mouse_control_samples <- fetch_control_samples('mouse')
human_control_samples <- fetch_control_samples('human')

# Function to export control sample IDs
export_control_sample_ids <- function(control_samples) {
  control_sample_ids <- control_samples$control_sample_ids
  write.csv(control_sample_ids, "control_sample_ids.csv", row.names = FALSE)
}

# Export control sample IDs
export_control_sample_ids(rat_control_samples)
export_control_sample_ids(mouse_control_samples)
export_control_sample_ids(human_control_samples)
