## Script to create the bundled pathway_res example dataset.
## Re-run this script (source it) whenever you want to regenerate the data object.
## Output: data/pathway_res.rda

library(readr)
library(dplyr)

# Read the full example file
full_data <- read_tsv(
  here::here("pathway-res-example.txt"),
  show_col_types = FALSE
)

# Retain only the columns used by reduceKappa_wrapper defaults, plus data_label
pathway_res <- full_data |>
  select(data_label, Geneset.ID, Description, Genes.Returned, P.value, FDR) |>
  # Keep the 50 most significant terms per data_label
  group_by(data_label) |>
  slice_min(FDR, n = 50, with_ties = FALSE) |>
  ungroup()

usethis::use_data(pathway_res, overwrite = TRUE)
