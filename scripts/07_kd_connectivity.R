#!/usr/bin/env Rscript
# Step 07 - LINCS gene-knockdown (KD) connectivity of predicted targets
# (Figure 4B).
#
# Two iLINCS queries against the LINCS consensus gene-knockdown library:
#   disease axis      : AMR signature (step 01)          vs KD library
#   atorvastatin axis : atorvastatin L1000 signature      vs KD library
#                       (LINCSCP_50693; PL21; 10 uM; 6 h; the top hit from step 02)
# Query signatures are restricted to the 100 most up- and 100 most
# down-regulated L1000 landmark genes, following the iLINCS connectivity
# framework (Pilarczyk et al., Nat Commun 2022). Similarity scores and P values
# are taken directly from iLINCS; when a gene has several KD signatures the
# median is used.
#
# Note on coverage: iLINCS returned only signatures with |similarity| >= 0.20
# for both queries (see data/cache/ilincs/), so a target whose knockdown has
# weaker connectivity on either axis has no score and is dropped.
#
# A target shows "opposing connectivity" when its disease and atorvastatin
# scores have opposite signs. r_disease < 0 with r_atorvastatin > 0 means the
# knockdown opposes AMR and resembles atorvastatin.
#
# Targets evaluated: the atorvastatin-only Venn region (63 targets) and the
# atorvastatin + rosuvastatin-only region (45 targets) from step 06.
#
# The cached queries in data/cache/ilincs/ are those reported in the
# manuscript. Set REFRESH_ILINCS=1 to re-query iLINCS.
#
# Outputs:
#   results/tables/kd_connectivity_atorvastatin_only_targets.csv
#   results/tables/kd_connectivity_atorvastatin_rosuvastatin_targets.csv

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(drugfindR)
  library(jsonlite)
})

n_per_side <- 100L
disease_cache <- file.path(PATHS$ilincs_cache, "kd_query_amr_signature.csv")
ator_cache <- file.path(PATHS$ilincs_cache, "kd_query_atorvastatin_signature.csv")

# Exact top-100 / bottom-100 logFC cutoffs on the L1000-matched signature
top_bottom_threshold <- function(dge, n = n_per_side) {
  sig <- prepareSignature(dge, geneColumn = "GeneSymbol", logfcColumn = "logFC", pvalColumn = "P.Value")
  values <- sort(sig$Value_LogDiffExp)
  n <- min(n, floor(length(values) / 2))
  c(values[n], values[length(values) - n + 1])
}

query_kd <- function(dge, source_name, similarity_threshold) {
  investigateSignature(
    expr = dge, outputLib = "KD",
    filterThreshold = top_bottom_threshold(dge),
    similarityThreshold = similarity_threshold, paired = TRUE,
    geneColumn = "GeneSymbol", logfcColumn = "logFC", pvalColumn = "P.Value",
    sourceName = source_name
  )
}

# ---- 1. Disease signature vs KD library ------------------------------------
if (file.exists(disease_cache) && !refresh_ilincs()) {
  message("Using cached disease-vs-KD query: ", disease_cache)
} else {
  dge <- read_csv(file.path(PATHS$intermediate, "dge_blood_amr_vs_dsa_negative.csv"), show_col_types = FALSE)
  write_csv(query_kd(dge, "Blood_DSA_rejection_vs_No_DSA", 0), disease_cache)
}

# ---- 2. Atorvastatin (LINCSCP_50693) signature vs KD library ---------------
if (file.exists(ator_cache) && !refresh_ilincs()) {
  message("Using cached atorvastatin-vs-KD query: ", ator_cache)
} else {
  x <- fromJSON(file.path(PATHS$raw_lincs, "LINCSCP_50693.json"))$data$signature
  ator_dge <- data.frame(
    GeneSymbol = toupper(trimws(x$Name_GeneSymbol)),
    logFC = as.numeric(x$Value_LogDiffExp),
    P.Value = as.numeric(x$Significance_pvalue)
  )
  ator_dge <- ator_dge[nzchar(ator_dge$GeneSymbol) & is.finite(ator_dge$logFC), ]
  write_csv(query_kd(ator_dge, "Atorvastatin_PL21", 0), ator_cache)
}

disease <- collapse_kd_by_gene(read_csv(disease_cache, show_col_types = FALSE)) |>
  rename(r_disease = Similarity, p_disease = pValue, n_disease = n_signatures)
ator <- collapse_kd_by_gene(read_csv(ator_cache, show_col_types = FALSE)) |>
  rename(r_ator = Similarity, p_ator = pValue, n_ator = n_signatures)

# ---- 3. Score each target set ----------------------------------------------
membership <- read_csv(file.path(PATHS$tables, "target_union_membership.csv"), show_col_types = FALSE)

score_targets <- function(region, out_name) {
  targets <- membership |> filter(Region == region)
  scored <- targets |>
    select(Target = TargetGene, AMRGene, AMRDirection) |>
    inner_join(disease, by = "Target") |>
    inner_join(ator, by = "Target") |>
    mutate(
      Relationship = if_else(sign(r_disease) != sign(r_ator), "Opposing", "Same direction"),
      OpposesAMR_ResemblesAtorvastatin = r_disease < 0 & r_ator > 0,
      BothNominalP05 = p_disease < 0.05 & p_ator < 0.05,
      Magnitude = abs(r_disease) + abs(r_ator)
    ) |>
    arrange(desc(Magnitude))
  write_csv(scored, file.path(PATHS$tables, out_name))
  message(sprintf(
    "%s: %d targets, %d with KD on both axes, %d opposing (%d oppose AMR + resemble atorvastatin), %d with P < 0.05 on both; AMR genes: %s",
    region, nrow(targets), nrow(scored), sum(scored$Relationship == "Opposing"),
    sum(scored$OpposesAMR_ResemblesAtorvastatin), sum(scored$BothNominalP05),
    paste(sprintf("%s (r_AMR %.2f, r_ator %.2f)", scored$Target[scored$AMRGene],
                  scored$r_disease[scored$AMRGene], scored$r_ator[scored$AMRGene]), collapse = ", ")
  ))
  scored
}

ator_only <- score_targets("Atorvastatin", "kd_connectivity_atorvastatin_only_targets.csv")
ator_rosu <- score_targets("Atorvastatin + Rosuvastatin", "kd_connectivity_atorvastatin_rosuvastatin_targets.csv")
message_done("KD connectivity tables written to ", PATHS$tables)
