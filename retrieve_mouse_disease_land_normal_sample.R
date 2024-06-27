library(opi)
library(stringr)
library(tidyverse)
library(openxlsx)

datastore <- opi()

# Function to randomly select 5 samples per tissue category
select_random_samples <- function(df, group_col, sample_size) {
  df %>%
    group_by({{group_col}}) %>%
    group_modify(~ {
      n_rows <- nrow(.x)
      if (n_rows <= sample_size) {
        .x
      } else {
        .x %>% slice_sample(n = sample_size)
      }
    }) %>%
    ungroup()
}

sql <- stringr::str_glue("
WITH numbered AS (
  SELECT
    sample_id, tissue, disease_state, spl.sample_index, ROW_NUMBER () OVER (PARTITION BY tissue ORDER BY rand ()) AS row_number
  FROM mouse_disease_b38.samples spl
  JOIN mouse_disease_b38.gene_fpkm fpkm on fpkm.sample_index = spl.sample_index
  WHERE spl.disease_state = 'normal control'  
)
SELECT sample_id,sample_index,tissue
FROM numbered 
WHERE row_number<2000
")
sample_list <- datastore$query(sql,dry_run=F)

# Apply the function to select 5 samples per tissue category
sample_list <- select_random_samples(sample_list, tissue, 3)

# View the result
print(sample_list)

unlisted_samples <- paste0(sample_list$sample_index,collapse=",")

sql <- stringr::str_glue("
SELECT sample_id,tissue,disease_state,fpkm,gene_name
FROM mouse_disease_b38.samples spl
JOIN mouse_disease_b38.gene_fpkm fpkm ON fpkm.sample_index = spl.sample_index
JOIN mouse_disease_b38.gene_annotation ann ON ann.gene_index = fpkm.gene_index
WHERE spl.sample_index IN ({unlisted_samples}) AND FPKM > 0
")
expression <- datastore$query(sql)
write.table(expression,file="mouse_normal_samples_from_mouse_diseaseland_20240607_N2.tsv",quote=F,sep="\t",row.names = F)
