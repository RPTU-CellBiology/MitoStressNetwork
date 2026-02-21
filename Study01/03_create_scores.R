library(here)
library(arrow)
library(dplyr)
library(tidyr)
library(msigdbr)
library(magrittr)
library(GSVA)
library(ggplot2)
library(ggridges)
library(forcats)
library(GGally)

data_dir <- here("Data", "CPTAC")

df_prot <- read_parquet(here(data_dir, "cptac_proteomics_filtered.parquet"))

# Check how many unmatched proteomics measurements are left
# TODO: Include missingness checks
unmatched <- df_prot |>
  add_count(Gene.ENSEMBL.ID, Sample.ID) |>
  filter(n != 2)

# Calculate Log2FC between matched samples for each protein
df_diff <- df_prot |>
  summarize(Protein.Abundance.Log2 = mean(Protein.Abundance.Log2, na.rm = TRUE),
            .by = c(Sample.ID, Sample.Type, Gene.ENSEMBL.ID)) |>
  pivot_wider(values_from = Protein.Abundance.Log2,
              names_from = Sample.Type,
              id_cols = c(Gene.ENSEMBL.ID, Sample.ID)) |>
  mutate(Log2FC = Tumor - Normal) |>
  drop_na(Log2FC)

# Enrichment Analysis (ssGSEA)
## Get gene sets
hallmark_gene_set <- msigdbr(species = "Homo sapiens", collection = "H")

gene_sets <- hallmark_gene_set |>
  distinct(gs_name, ensembl_gene) |>
  filter(!is.na(ensembl_gene)) |>
  group_by(gs_name) |>
  summarise(genes = list(unique(ensembl_gene)), .groups = "drop") %$%
  setNames(genes, gs_name)

## Transform data to matrix
mat_diff <- df_diff |>
  select(Sample.ID, Gene.ENSEMBL.ID, Log2FC) |>
  pivot_wider(values_from = Log2FC, names_from = Sample.ID, id_cols = Gene.ENSEMBL.ID) |>
  tibble::column_to_rownames("Gene.ENSEMBL.ID") |>
  na.omit() |>
  as.matrix()

## Perform ssGSEA
gseaPar <- GSVA::ssgseaParam(mat_diff, gene_sets)
gsea_result <- GSVA::gsva(gseaPar)

df_gsea <- gsea_result |>
  as.data.frame() |>
  tibble::rownames_to_column("GeneSet") |>
  pivot_longer(-GeneSet, names_to = "Sample.ID", values_to = "ssGSEA.Score") |>
  write_parquet(here(data_dir, "ssGSEA.parquet"))

## Plot enrichment score distributions
plot_enrichment_scores <- function (df) {
  colors <- rev(RColorBrewer::brewer.pal(11, "RdBu"))

  df |>
    ggplot() +
    aes(x = ssGSEA.Score, y = fct_reorder(GeneSet, -ssGSEA.Score), fill = after_stat(x)) +
    geom_vline(xintercept = 0) +
    geom_density_ridges_gradient() +
    scale_fill_gradientn(colors = colors, space = "Lab",
                         limits = c(-0.5, 0.5), oob = scales::squish, guide = NULL) +
    labs(x = "ssGSEA Score", y = "Gene Set") +
    theme_bw()
}

plot_enrichment_scores(df_gsea)

# Calculate Proxy Scores
df_scores <- df_gsea |>
  mutate(ssGSEA.Score.Z = scale(ssGSEA.Score)[,1], .by = GeneSet) |>
  summarize(
    Score.Prolif = mean(ssGSEA.Score.Z[GeneSet %in% c(
      "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"
    )]),
    Score.Hypoxia = mean(ssGSEA.Score.Z[GeneSet == "HALLMARK_HYPOXIA"]),
    Score.Glycolysis = mean(ssGSEA.Score.Z[GeneSet == "HALLMARK_GLYCOLYSIS"]),
    Score.MitoFitness = mean(ssGSEA.Score.Z[GeneSet == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]),
    .by = Sample.ID) |>
  write_parquet(here(data_dir, "sample_scores.parquet"))

## Correlate Scores
df_scores |>
  select(starts_with("Score.")) |>
  ggpairs()

## Plot GSEA scores of samples with low mito fitness
df_gsea |>
  semi_join(df_scores |> filter(Score.MitoFitness < -0.2), by = "Sample.ID") |>
  plot_enrichment_scores()

# TODO: Consider adding MYC targets to proliferation scores (also cancer stage and other clinical indicators)
# TODO: Add further gene sets to mito fitness scores (GO_MITOCHONDRIAL_TRANSLATION, GO_MITOCHONDRIAL_RIBOSOME, REACTOME_RESPIRATORY_ELECTRON_TRANSPORT)
# TODO: Consider imputation instead of na.omit
# TODO: Check coverage of IDs between gene sets and dataset
# TODO: Consider combining ssGSEA scores via PCA
# TODO: Scale scores
# TODO: Keep GO_RESPONSE_TO_MITOCHONDRIAL_STRESS as validation score for stress network
# TODO: Validate proliferation score (correlation: MKI67 (Ki-67) (if present), PCNA, MCM2–MCM7, TOP2A, CDK1, CCNB1/2, CDC20, UBE2C)
# TODO: Validate MitoFitness score (correlation: VDAC1, TOMM20, TOMM40, NDUFS, SDHB, UQCRC1/2, COX4I1, ATP5F1A/B)
# TODO: Stratify validation across cancer types; consider bootstrapping
