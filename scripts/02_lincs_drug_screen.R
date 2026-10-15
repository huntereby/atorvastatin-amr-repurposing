#!/usr/bin/env Rscript
# Step 02 - Signature-based drug-repurposing screen (LINCS L1000 / iLINCS).
#
# The whole-blood AMR signature from step 01 is reduced to the L1000 landmark
# genes, the top and bottom 10% by log2 fold change are kept, and the up/down
# lists are queried (paired) against the iLINCS chemical-perturbagen (CP)
# library with drugfindR. Negative similarity = discordant = the compound's
# signature opposes the AMR signature.
#
# Discordant compounds are ranked from most to least negative similarity.
# Drug ranks refer to position among ALL discordant compounds. FDA status is
# taken from Drugs@FDA (currently marketed prescription or OTC products) and
# mechanisms of action for the 100 most discordant compounds from the Broad
# Drug Repurposing Hub (release 2025-08-18).
#
# The iLINCS query is read from data/cache/ilincs/cp_query_amr_signature.csv
# (the result reported in the manuscript). Set REFRESH_ILINCS=1 to re-query.
#
# Outputs:
#   results/intermediate/lincs_cp_results.csv          all compounds (case-collapsed)
#   results/tables/lincs_discordant_fda_approved.csv   FDA-approved discordant drugs
#   results/tables/lincs_top10_fda_approved.csv        Figure 1 data
#   results/tables/lincs_top100_moa_annotations.csv
#   results/tables/lincs_top100_moa_counts.csv

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(drugfindR)
  library(stringr)
  library(tidyr)
})

filter_prop <- 0.10
similarity_threshold <- 0.20
source_name <- "Blood_DSA_rejection_vs_No_DSA"

dge_file <- file.path(PATHS$intermediate, "dge_blood_amr_vs_dsa_negative.csv")
cache_file <- file.path(PATHS$ilincs_cache, "cp_query_amr_signature.csv")

# ---------------------------------------------------------------------------
# 1. iLINCS chemical-perturbagen query
# ---------------------------------------------------------------------------
if (file.exists(cache_file) && !refresh_ilincs()) {
  message("Using cached iLINCS CP query: ", cache_file)
  cp <- read_csv(cache_file, show_col_types = FALSE)
} else {
  dge <- read_csv(dge_file, show_col_types = FALSE)
  cp <- investigateSignature(
    expr = dge,
    outputLib = "CP",
    filterProp = filter_prop,
    similarityThreshold = similarity_threshold,
    paired = TRUE,
    geneColumn = "GeneSymbol",
    logfcColumn = "logFC",
    pvalColumn = "P.Value",
    sourceName = source_name
  )
  write_csv(cp, cache_file)
}

# iLINCS has a few treatment labels that differ only by case (e.g.
# "atorvastatin" / "Atorvastatin"). Keep the signature with the largest
# absolute similarity for each, matching drugfindR's per-treatment rule.
cp <- cp |>
  mutate(Target = str_squish(Target), TargetKey = str_to_lower(Target)) |>
  group_by(TargetKey) |>
  slice_max(abs(Similarity), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(Similarity)
write_csv(cp, file.path(PATHS$intermediate, "lincs_cp_results.csv"))

discordant <- cp |>
  filter(Similarity < 0) |>
  arrange(Similarity) |>
  mutate(OverallDiscordantRank = row_number())

message(sprintf("%d compounds returned; %d discordant. Top hit: %s (%.2f; %s, %s, %s, %s)",
                nrow(cp), nrow(discordant), discordant$Target[1], discordant$Similarity[1],
                discordant$TargetSignature[1], discordant$TargetCellLine[1],
                discordant$TargetConcentration[1], discordant$TargetTime[1]))

# ---------------------------------------------------------------------------
# 2. FDA-approved, currently marketed (Drugs@FDA)
# ---------------------------------------------------------------------------
normalize_key <- function(x) str_replace_all(str_to_lower(x), "[^a-z0-9]", "")

fda_dir <- file.path(PATHS$external_cache, "drugsatfda")
if (!file.exists(file.path(fda_dir, "Products.txt"))) {
  dir.create(fda_dir, recursive = TRUE, showWarnings = FALSE)
  fda_zip <- file.path(fda_dir, "drugsatfda.zip")
  download.file("https://www.fda.gov/media/89850/download?attachment=", fda_zip, mode = "wb", quiet = TRUE)
  unzip(fda_zip, exdir = fda_dir)
}
read_fda <- function(f) {
  read.delim(file.path(fda_dir, f), sep = "\t", quote = "", fill = TRUE,
             check.names = FALSE, stringsAsFactors = FALSE)
}

fda_products <- read_fda("Products.txt") |>
  inner_join(read_fda("MarketingStatus.txt"), by = c("ApplNo", "ProductNo")) |>
  left_join(read_fda("MarketingStatus_Lookup.txt"), by = "MarketingStatusID") |>
  filter(MarketingStatusID %in% c(1, 2)) |> # 1 = prescription, 2 = OTC
  mutate(DrugNameKey = normalize_key(DrugName), ActiveIngredientKey = normalize_key(ActiveIngredient))

match_fda_product <- function(target) {
  key <- normalize_key(target)
  hits <- fda_products |>
    mutate(MatchPriority = case_when(
      ActiveIngredientKey == key ~ 1L,
      DrugNameKey == key ~ 2L,
      str_starts(ActiveIngredientKey, fixed(key)) ~ 3L,
      str_starts(DrugNameKey, fixed(key)) ~ 4L,
      TRUE ~ 99L
    )) |>
    filter(MatchPriority < 99L) |>
    arrange(MatchPriority, desc(ReferenceDrug), ApplNo, ProductNo)
  if (!nrow(hits)) return(NULL)
  hits |>
    slice(1) |>
    transmute(
      Target = target,
      FDAProduct = str_to_title(DrugName),
      FDAActiveIngredient = str_to_title(ActiveIngredient),
      FDAApplication = ApplNo,
      FDAMarketingStatus = MarketingStatusDescription,
      FDAMatchType = c("Exact active ingredient", "Exact product/brand name",
                       "Active ingredient salt/form", "Product/brand salt/form")[MatchPriority]
    )
}

fda_discordant <- discordant |>
  inner_join(bind_rows(lapply(discordant$Target, match_fda_product)), by = "Target") |>
  arrange(OverallDiscordantRank)
write_csv(fda_discordant, file.path(PATHS$tables, "lincs_discordant_fda_approved.csv"))

top_fda <- fda_discordant |>
  slice_head(n = 10) |>
  select(OverallDiscordantRank, Target, Similarity, TargetSignature, TargetCellLine,
         TargetConcentration, TargetTime, pValue, FDAProduct, FDAActiveIngredient, FDAMarketingStatus)
write_csv(top_fda, file.path(PATHS$tables, "lincs_top10_fda_approved.csv"))

# ---------------------------------------------------------------------------
# 3. Mechanisms of action for the 100 most discordant compounds
# ---------------------------------------------------------------------------
broad_file <- file.path(PATHS$external_cache, "repo-drug-annotation-20250818.txt")
if (!file.exists(broad_file)) {
  download.file("https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20250818.txt",
                broad_file, mode = "wb", quiet = TRUE)
}
broad <- read.delim(broad_file, sep = "\t", quote = "\"", skip = 9, fill = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE) |>
  mutate(AnnotationKey = normalize_key(pert_iname)) |>
  group_by(AnnotationKey) |>
  summarise(
    BroadName = first(pert_iname),
    ClinicalPhase = paste(sort(unique(clinical_phase[nzchar(clinical_phase)])), collapse = " | "),
    MOA = paste(sort(unique(moa[nzchar(moa)])), collapse = " | "),
    .groups = "drop"
  )

# iLINCS names that differ from the Repurposing Hub name for the same compound
annotation_aliases <- c(
  "amidate" = "etomidate",
  "9fluoro111721trihydroxy16methylpregna14diene320dione" = "dexamethasone",
  "quinidinegluconate" = "quinidine"
)

top100 <- discordant |>
  slice_head(n = 100) |>
  mutate(
    AnnotationKey = normalize_key(Target),
    AnnotationKey = if_else(AnnotationKey %in% names(annotation_aliases),
                            unname(annotation_aliases[AnnotationKey]), AnnotationKey)
  ) |>
  left_join(broad, by = "AnnotationKey") |>
  mutate(HasCuratedMOA = !is.na(MOA) & nzchar(MOA))
write_csv(
  top100 |> select(OverallDiscordantRank, Target, Similarity, TargetSignature, TargetCellLine,
                   BroadName, ClinicalPhase, MOA, HasCuratedMOA),
  file.path(PATHS$tables, "lincs_top100_moa_annotations.csv")
)

pluralize <- function(x) {
  if_else(str_detect(x, "(inhibitor|antagonist|agonist|activator|modulator|blocker|agent|scavenger)$"),
          paste0(x, "s"), x)
}
moa_counts <- top100 |>
  filter(HasCuratedMOA) |>
  separate_rows(MOA, sep = "\\s*\\|\\s*") |>
  mutate(
    MOALabel = pluralize(str_to_sentence(str_squish(MOA))),
    MOALabel = str_replace_all(MOALabel, c(
      "\\bParp\\b" = "PARP", "\\bPi3k\\b" = "PI3K", "\\bBcl\\b" = "BCL",
      "\\bP38 mapk\\b" = "p38 MAPK", "\\bHmgcr\\b" = "HMGCR", "\\bFlt3\\b" = "FLT3",
      "\\bKit\\b" = "KIT", "\\bEgfr\\b" = "EGFR", "\\bVegfr\\b" = "VEGFR",
      "\\bAkt\\b" = "AKT", "\\bAtr\\b" = "ATR", "\\bNfkb\\b" = "NF-kB", "\\bStat\\b" = "STAT"
    ))
  ) |>
  distinct(OverallDiscordantRank, Target, MOALabel, .keep_all = TRUE) |>
  group_by(MOALabel) |>
  summarise(
    DrugCount = n_distinct(Target),
    BestDiscordantRank = min(OverallDiscordantRank),
    Drugs = paste(Target[order(OverallDiscordantRank)], collapse = "; "),
    .groups = "drop"
  ) |>
  arrange(desc(DrugCount), BestDiscordantRank, MOALabel)
write_csv(moa_counts, file.path(PATHS$tables, "lincs_top100_moa_counts.csv"))

message_done(sprintf("%d FDA-approved discordant drugs; %d of top 100 have a curated MOA",
                     nrow(fda_discordant), sum(top100$HasCuratedMOA)))
