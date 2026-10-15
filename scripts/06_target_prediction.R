#!/usr/bin/env Rscript
# Step 06 - Predicted protein targets of four statins and overlap with the
# AMR gene set (Figure 4A).
#
# Inputs (data/raw/target_prediction/), one file per statin:
#   swisstargetprediction/<statin>_swisstargetprediction.csv  top 100 human targets
#   superpred/<statin>_superpred.csv                          all returned targets
#     (SuperPred web server output; every retained prediction has
#      probability >= 50% and model accuracy >= 70%)
#
# UniProt accessions are mapped to HGNC symbols with org.Hs.eg.db. For
# SwissTargetPrediction, the "Common name" is used when no mapping exists.
# Each statin's target set is the union of the two methods. Targets are then
# assigned to the 15 mutually exclusive regions of the four-set Venn diagram
# and intersected with the AMR gene set from step 01.
#
# Outputs:
#   results/tables/target_predictions_per_statin.csv   per statin x target, with method evidence
#   results/tables/target_union_membership.csv         431 targets x 4 statins (0/1)
#   results/tables/target_venn_regions_amr.csv         15 regions: size, AMR hits, genes
#   results/tables/method_overlap_per_statin.csv       Swiss vs SuperPred counts per statin

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(dplyr)
  library(org.Hs.eg.db)
  library(purrr)
  library(stringr)
  library(tidyr)
})

drugs <- STATINS_CLINICAL
swiss_dir <- file.path(PATHS$raw_targets, "swisstargetprediction")
superpred_dir <- file.path(PATHS$raw_targets, "superpred")

clean_gene <- function(x) {
  x |> str_trim() |> str_to_upper() |> na_if("") |> na_if("NA") |> na_if("NOT AVAILABLE")
}

uniprot_to_symbol <- function(ids) {
  ids <- ids |> str_trim() |> na_if("") |> na_if("Not Available")
  keys <- unique(na.omit(ids))
  if (!length(keys)) return(rep(NA_character_, length(ids)))
  mapped <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = keys, keytype = "UNIPROT",
                                  column = "SYMBOL", multiVals = "first")
  unname(mapped[ids])
}

read_swiss <- function(drug) {
  read_csv(file.path(swiss_dir, paste0(tolower(drug), "_swisstargetprediction.csv")),
           show_col_types = FALSE) |>
    mutate(
      Drug = drug,
      UniProtID = `Uniprot ID`,
      TargetGene = coalesce(clean_gene(uniprot_to_symbol(`Uniprot ID`)), clean_gene(`Common name`)),
      TargetName = Target,
      SwissProbability = `Probability*`
    ) |>
    filter(!is.na(TargetGene)) |>
    distinct(Drug, TargetGene, .keep_all = TRUE) |>
    select(Drug, TargetGene, UniProtID, TargetName, SwissProbability)
}

read_superpred <- function(drug) {
  read_csv(file.path(superpred_dir, paste0(tolower(drug), "_superpred.csv")),
           show_col_types = FALSE) |>
    mutate(
      Drug = drug,
      UniProtID = `UniProt ID`,
      TargetGene = clean_gene(uniprot_to_symbol(`UniProt ID`)),
      TargetName = `Target Name`,
      SuperPredProbability = parse_number(Probability) / 100,
      SuperPredModelAccuracy = parse_number(`Model accuracy`) / 100
    ) |>
    filter(!is.na(TargetGene)) |>
    distinct(Drug, TargetGene, .keep_all = TRUE) |>
    select(Drug, TargetGene, UniProtID, TargetName, SuperPredProbability, SuperPredModelAccuracy)
}

swiss <- map_dfr(drugs, read_swiss)
superpred <- map_dfr(drugs, read_superpred)
stopifnot(all(superpred$SuperPredProbability >= 0.5), all(superpred$SuperPredModelAccuracy >= 0.7))

per_statin <- full_join(
  swiss |> rename(SwissUniProtID = UniProtID, SwissTargetName = TargetName),
  superpred |> rename(SuperPredUniProtID = UniProtID, SuperPredTargetName = TargetName),
  by = c("Drug", "TargetGene")
) |>
  mutate(Evidence = case_when(
    !is.na(SwissProbability) & !is.na(SuperPredProbability) ~ "SwissTargetPrediction + SuperPred",
    !is.na(SwissProbability) ~ "SwissTargetPrediction only",
    TRUE ~ "SuperPred only"
  )) |>
  arrange(Drug, TargetGene)
write_csv(per_statin, file.path(PATHS$tables, "target_predictions_per_statin.csv"))

method_overlap <- per_statin |>
  group_by(Drug) |>
  summarise(
    SwissTargets = sum(!is.na(SwissProbability)),
    SuperPredTargets = sum(!is.na(SuperPredProbability)),
    BothMethods = sum(Evidence == "SwissTargetPrediction + SuperPred"),
    UnionTargets = n(),
    .groups = "drop"
  ) |>
  arrange(match(Drug, drugs))
write_csv(method_overlap, file.path(PATHS$tables, "method_overlap_per_statin.csv"))

# ---- Four-statin membership and Venn regions --------------------------------
amr <- read_csv(file.path(PATHS$intermediate, "amr_gene_set.csv"), show_col_types = FALSE) |>
  arrange(adj.P.Val)

membership <- per_statin |>
  distinct(TargetGene, Drug) |>
  mutate(In = 1L) |>
  pivot_wider(names_from = Drug, values_from = In, values_fill = 0L) |>
  select(TargetGene, all_of(drugs)) |>
  mutate(
    StatinCount = rowSums(pick(all_of(drugs))),
    Region = apply(pick(all_of(drugs)), 1, function(r) paste(drugs[r == 1], collapse = " + ")),
    AMRGene = TargetGene %in% amr$GeneSymbol
  ) |>
  left_join(amr |> select(TargetGene = GeneSymbol, AMRDirection = Direction,
                          AMR_logFC = logFC, AMR_FDR = adj.P.Val),
            by = "TargetGene") |>
  arrange(desc(StatinCount), TargetGene)
write_csv(membership, file.path(PATHS$tables, "target_union_membership.csv"))

all_regions <- unlist(lapply(seq_along(drugs), function(k) combn(drugs, k, simplify = FALSE)),
                      recursive = FALSE)
regions <- map_dfr(all_regions, function(in_drugs) {
  label <- paste(in_drugs, collapse = " + ")
  genes <- membership |> filter(Region == label)
  hits <- genes |> filter(AMRGene) |> arrange(AMR_FDR) # most significant first
  tibble(
    Region = label,
    nStatins = length(in_drugs),
    Targets = nrow(genes),
    AMRGenes = nrow(hits),
    AMRGeneList = paste(hits$TargetGene, collapse = ";")
  )
})
write_csv(regions, file.path(PATHS$tables, "target_venn_regions_amr.csv"))

message(paste(sprintf("%-13s %d targets", drugs,
                      sapply(drugs, function(d) sum(membership[[d]]))), collapse = "\n"))
message_done(sprintf(
  "%d unique targets; %d shared by all four; %d (%.1f%%) in the AMR gene set; atorvastatin+rosuvastatin-only region: %s",
  nrow(membership), sum(membership$StatinCount == 4), sum(membership$AMRGene),
  100 * mean(membership$AMRGene),
  with(filter(regions, Region == "Atorvastatin + Rosuvastatin"), sprintf("%d targets, AMR genes %s", Targets, AMRGeneList))
))
