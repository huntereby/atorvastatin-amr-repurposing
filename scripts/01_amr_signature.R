#!/usr/bin/env Rscript
# Step 01 - Whole-blood AMR transcriptional signature (GSE50084).
#
# Downloads GSE50084 (Affymetrix Human Gene 1.0 ST, GPL6244) from GEO,
# collapses probes to genes (highest mean expression probe per gene), and fits
# a no-intercept limma model comparing whole-blood samples from DSA-positive
# recipients with biopsy-proven AMR (n = 28) against DSA-negative recipients
# (n = 12). P values are Benjamini-Hochberg adjusted.
#
# The AMR gene set used for target overlap (Methods; Figure S1) is the top
# ceiling(10%) of genes in each direction ranked by log2 fold change that also
# reach FDR < 0.05.
#
# Outputs (results/intermediate/):
#   gse50084_sample_metadata.csv
#   dge_blood_amr_vs_dsa_negative.csv     full limma table (input to steps 02, 08)
#   amr_top10pct_up.csv / amr_top10pct_down.csv
#   amr_gene_set.csv                      top 10% AND FDR < 0.05 (2,046 genes)

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(Biobase)
  library(dplyr)
  library(GEOquery)
  library(limma)
  library(stringr)
  library(tibble)
})

case_group <- "DSA_rejection"
control_group <- "No_DSA"
top_prop <- 0.10
fdr_cutoff <- 0.05

message("Downloading GSE50084 series matrix...")
gse <- getGEO("GSE50084", GSEMatrix = TRUE, getGPL = FALSE)[[1]]
expr_probe <- exprs(gse)

pheno <- pData(gse) |>
  rownames_to_column("sample_id") |>
  mutate(
    sample_title = title,
    tissue_group = case_when(
      str_detect(sample_title, "^Biopsy_") ~ "Biopsy",
      str_detect(sample_title, "^Blood_") ~ "Blood",
      TRUE ~ NA_character_
    ),
    disease_group = case_when(
      str_detect(sample_title, "_DSA_rejection_") ~ "DSA_rejection",
      str_detect(sample_title, "_DSA_normal_") ~ "DSA_normal",
      str_detect(sample_title, "_No_DSA_") ~ "No_DSA",
      TRUE ~ NA_character_
    )
  )
write_csv(
  pheno |> select(sample_id, geo_accession, sample_title, tissue_group, disease_group),
  file.path(PATHS$intermediate, "gse50084_sample_metadata.csv")
)

message("Downloading GPL6244 annotation...")
gpl <- getGEO("GPL6244")
annot <- Table(gpl) |>
  transmute(
    ID = as.character(ID),
    GeneSymbol = str_match(gene_assignment, " // ([^/]+?) // ")[, 2]
  ) |>
  mutate(
    GeneSymbol = str_trim(GeneSymbol),
    GeneSymbol = na_if(GeneSymbol, "---"),
    GeneSymbol = na_if(GeneSymbol, "")
  )

# One row per gene: keep the probe set with the highest mean expression.
sample_cols <- colnames(expr_probe)
expr_gene <- as.data.frame(expr_probe) |>
  rownames_to_column("ID") |>
  left_join(annot, by = "ID") |>
  filter(!is.na(GeneSymbol)) |>
  mutate(.mean_expr = rowMeans(across(all_of(sample_cols)), na.rm = TRUE)) |>
  group_by(GeneSymbol) |>
  arrange(desc(.mean_expr), .by_group = TRUE) |>
  slice(1) |>
  ungroup() |>
  select(GeneSymbol, all_of(sample_cols)) |>
  column_to_rownames("GeneSymbol") |>
  as.matrix()
storage.mode(expr_gene) <- "numeric"

# ---- limma: whole blood, AMR (DSA rejection) vs DSA-negative ---------------
keep <- pheno |>
  filter(tissue_group == "Blood", disease_group %in% c(case_group, control_group))
stopifnot(sum(keep$disease_group == case_group) == 28, sum(keep$disease_group == control_group) == 12)

x <- expr_gene[, keep$geo_accession, drop = FALSE]
group <- factor(keep$disease_group, levels = c(control_group, case_group))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
contrast <- makeContrasts(contrasts = paste(case_group, control_group, sep = "-"), levels = design)

fit <- eBayes(contrasts.fit(lmFit(x, design), contrast))
dge <- topTable(fit, coef = 1, number = Inf, sort.by = "P") |>
  rownames_to_column("GeneSymbol")
write_csv(dge, file.path(PATHS$intermediate, "dge_blood_amr_vs_dsa_negative.csv"))

# ---- AMR gene set -----------------------------------------------------------
n_top <- ceiling(nrow(dge) * top_prop)
top_up <- dge |> arrange(desc(logFC)) |> slice_head(n = n_top)
top_down <- dge |> arrange(logFC) |> slice_head(n = n_top)
write_csv(top_up, file.path(PATHS$intermediate, "amr_top10pct_up.csv"))
write_csv(top_down, file.path(PATHS$intermediate, "amr_top10pct_down.csv"))

amr_gene_set <- bind_rows(
  top_up |> mutate(Direction = "Up"),
  top_down |> mutate(Direction = "Down")
) |>
  filter(adj.P.Val < fdr_cutoff) |>
  mutate(GeneSymbol = toupper(trimws(GeneSymbol))) |>
  distinct(GeneSymbol, .keep_all = TRUE) |>
  select(GeneSymbol, Direction, logFC, P.Value, adj.P.Val)
write_csv(amr_gene_set, file.path(PATHS$intermediate, "amr_gene_set.csv"))

message_done(sprintf(
  "%d genes tested; %d at FDR < 0.05; AMR gene set = %d (%d up, %d down); fold-change cutoffs %.2f / +%.2f",
  nrow(dge), sum(dge$adj.P.Val < fdr_cutoff), nrow(amr_gene_set),
  sum(amr_gene_set$Direction == "Up"), sum(amr_gene_set$Direction == "Down"),
  max(top_down$logFC), min(top_up$logFC)
))
