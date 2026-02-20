library(here)
library(arrow)
library(dplyr)
library(tidyr)

data_dir <- here("Data", "CPTAC")

df_prot <- read_parquet(here(data_dir, "cptac_proteomics_filtered.parquet"))

# Check how many unmatched proteomics measurements are left
# TODO: Include missingness checks
df_prot |>
  add_count(idx, Sample.ID) |>
  filter(n < 2) |>
  nrow()

# Calculate Log2FC between matched samples for each protein
df_diff <- df_prot |>
  pivot_wider(values_from = Protein.Abundance.Log2,
              names_from = Sample.Type,
              id_cols = c(idx, Sample.ID)) |>
  mutate(Log2FC = Tumor - Normal) |>
  drop_na(Log2FC)

# TODO: Perform GSEA
# TODO: Calculate proxy scores per sample for OXPHOS, Prolfieration, and Hypoxia
# TODO: Correlate proxy scores
# TODO: Remove samples that don't have significantly enriched OXPHOS signals
# TODO: Calculate gene-wise pairwise partial correlation (covariates: purity, ploidy, cancer type, proliferation, gene copy number difference)
#
