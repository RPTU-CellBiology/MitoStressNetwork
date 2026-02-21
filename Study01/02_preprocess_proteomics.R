library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(vroom)
library(here)
library(arrow)

here::i_am(".Rprofile")

data_dir <- here("Data", "CPTAC")
prot_data_dir <- here(data_dir, "Proteomics", "Proteome_BCM_GENCODE_v34_harmonized_v1")
meta_data_dir <- here(prot_data_dir, "README")

# === Load Datasets ===
## Proteomics
read_prot <- function(filename, directory = prot_data_dir) {
  vroom::vroom(here(directory, filename), show_col_types = FALSE) |>
    pivot_longer(-idx, names_to = "Sample.ID", values_to = "Protein.Abundance.Log2",
                 values_transform = as.numeric)
}

df_prot <- list.files(prot_data_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, "_proteomics_.*\\.txt$")) |>
  tibble(Filename = _) |>
  mutate(Sample.Type = str_split_i(str_remove(Filename, "\\.txt$"), "_", 9)) |>
  mutate(data = map(Filename, read_prot)) |>
  unnest(data) |>
  mutate(Gene.ENSEMBL.ID = sub("(\\.\\d+)$", "", idx)) |>
  select(Sample.ID, Sample.Type, Gene.ENSEMBL.ID, Protein.Abundance.Log2)

## Metadata
read_meta <- function(filename, directory = meta_data_dir) {
  vroom::vroom(here(directory, filename), comment = "data_type", show_col_types = FALSE)
}

df_meta <- list.files(meta_data_dir, recursive = FALSE, full.names = FALSE) |>
  keep(~ str_detect(.x, "_meta\\.txt$")) |>
  tibble(Filename = _) |>
  mutate(Sample.CancerType = str_split_i(sub('\\.txt$', '', Filename), "_", 1)) |>
  mutate(data = map(Filename, read_meta)) |>
  unnest(data) |>
  select(-Filename) |>
  rename(Sample.ID = idx)

# === Filter Proteomics ===

df_meta_filtered <- df_meta |>
  mutate(Sample.TumorPurity = pmin(WES_purity, WGS_purity, na.rm = TRUE)) |>
  filter(Tumor == "Yes" & Normal == "Yes") |>
  filter(Sample.TumorPurity > 0.2) |>
  write_parquet(here(data_dir, 'cptac_meta_filtered.parquet'))

df_prot_filtered <- df_prot |>
  semi_join(y = df_meta_filtered, by = "Sample.ID") |>
  write_parquet(here(data_dir, 'cptac_proteomics_filtered.parquet'))
