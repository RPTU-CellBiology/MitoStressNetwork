library(here)
library(arrow)
library(dplyr)
library(limma)
library(WGCNA)

here::i_am(".Rprofile")

data_dir <- here("Data", "CPTAC")

df_diff <- read_parquet(here(data_dir, "cptac_diff.parquet"))
df_scores <- read_parquet(here(data_dir, "sample_scores.parquet"))
df_gsea <- read_parquet(here(data_dir, "ssGSEA.parquet"))
df_meta <- read_parquet(here(data_dir, 'cptac_meta_filtered.parquet'))

# Remove samples with high mito fitness
df_scores_filtered <- df_scores |>
  filter(Score.MitoFitness < -1)

# Combine confounding variables (purity, ploidy, cancer type, p53 status, proliferation, gene copy number difference)
df_conf <- df_meta |>
  mutate(Sample.Ploidy = median(WGS_ploidy, WES_ploidy, na.rm = TRUE),
         # Median imputation of ploidy
         Sample.Ploidy = if_else(is.na(Sample.Ploidy), median(Sample.Ploidy, na.rm = TRUE), Sample.Ploidy),
         .by = Sample.ID) |>
  mutate(Sample.TumorPurity = Sample.TumorPurity,
         # Sample.TP53mut = first(if_else(is.na(TP53_mutation), "Unknown", as.character(as.logical(TP53_mutation)))),
         Score.TP53 = scale(HALLMARK_P53_PATHWAY)[, 1]) |> # Transcriptomics-based ssGSEAc TP53 pathway scores as indicator for p53 activity
  select(Sample.ID, Sample.TumorPurity, Sample.CancerType, Sample.Ploidy, Score.TP53) |>
  drop_na() |>
  inner_join(y = df_scores_filtered, by = "Sample.ID") |>
  mutate(Sample.CancerType = factor(Sample.CancerType))

# Build linear model to residualize fold-changes
mat_fc <- df_diff |>
  select(Gene.ENSEMBL.ID, Sample.ID, Log2FC) |>
  semi_join(y = df_conf, by = "Sample.ID") |>
  pivot_wider(names_from = Sample.ID, values_from = Log2FC, id_cols = Gene.ENSEMBL.ID) |>
  tibble::column_to_rownames("Gene.ENSEMBL.ID") |>
  na.omit() |>
  as.matrix()

design <- model.matrix(
  ~ 1 + Sample.Ploidy + Sample.TumorPurity + Sample.CancerType + Score.TP53 + Score.Prolif,
  data = df_conf
)

fit <- lmFit(mat_fc, design)
fitted <- fit$coefficients %*% t(design)   # genes x samples
mat_resid <- mat_fc - fitted

# Calculate gene-wise co-expression network using residuals of linear model
options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()

datExpr <- t(mat_resid)  # samples x genes

## Remove bad genes/samples
gsg <- goodSamplesGenes(datExpr, verbose = 3)
datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]

## Choose soft-threshold
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5)

## Build network
net <- blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  networkType = "signed",
  TOMType = "signed",
  corType = "bicor",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

module_labels <- net$colors
module_eigengenes <- net$MEs

trait <- df_conf$Score.MitoFitness
fivenum(sapply(module_eigengenes, \(x) cor(x, trait)))

df_modules <- module_labels |>
  tibble::enframe()

# TODO: Consider whether to include hypoxia/glycolysis as further covariates
# TODO: Consider tumor stage/histologic grade/PFS (prolif marker, covariate, phenotypic analysis)
# TODO: Check if sample IDs match for each step
# TODO: Map to HGNC IDs for better interpretation
