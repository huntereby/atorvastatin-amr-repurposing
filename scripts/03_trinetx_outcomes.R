#!/usr/bin/env Rscript
# Step 03 - One-year kidney transplant outcomes in TriNetX.
#
# Reads the TriNetX "Compare Outcomes" exports in data/raw/trinetx/ (1:1
# propensity-score-matched cohorts; outcome window day 1-365 after renal
# allotransplantation, CPT 1008109):
#   Outcome 1 = kidney transplant rejection (ICD-10-CM T86.11)
#   Outcome 2 = kidney transplant failure   (ICD-10-CM T86.12)
# (TriNetX labels both "Unnamed Outcome" in the export metadata; the outcome
# numbering above follows the analysis set-up.)
#
# Risk ratios, risk differences and their 95% CIs are taken directly from
# TriNetX. NNT = 1 / absolute risk reduction, with its CI obtained by
# inverting the risk-difference CI limits. NNTs are stored unrounded and
# reported rounded up to the next whole patient.
#
# Outputs:
#   results/tables/trinetx_each_statin_vs_no_statin.csv   (Figure 2A)
#   results/tables/trinetx_atorvastatin_vs_other_statins.csv (Figure 2B)

source(file.path("R", "utils.R"))
suppressPackageStartupMessages(library(dplyr))

direction <- function(df) {
  case_when(
    df$RR_Upper < 1 ~ "Lower risk",
    df$RR_Lower > 1 ~ "Higher risk",
    TRUE ~ "Not significant"
  )
}

read_study <- function(study) {
  study_dir <- file.path(PATHS$raw_trinetx, study)
  bind_rows(read_trinetx_outcome(study_dir, 1), read_trinetx_outcome(study_dir, 2))
}

# ---- Primary analysis: each statin vs statin-naive -------------------------
each_statin <- bind_rows(lapply(STATINS_CLINICAL, function(statin) {
  read_study(paste0(tolower(statin), "_vs_no_statin")) |>
    mutate(Statin = statin, .before = 1)
}))
each_statin$Direction <- direction(each_statin)
each_statin <- add_nnt(each_statin)
write_csv(each_statin, file.path(PATHS$tables, "trinetx_each_statin_vs_no_statin.csv"))

# ---- Active comparator: atorvastatin vs other statins ----------------------
active <- read_study("atorvastatin_vs_other_statins") |>
  mutate(Comparison = "Atorvastatin vs other statins", .before = 1)
active$Direction <- direction(active)
active <- add_nnt(active)
write_csv(active, file.path(PATHS$tables, "trinetx_atorvastatin_vs_other_statins.csv"))

fmt <- function(df, label) {
  sprintf("%-14s %-30s n=%5d/%5d  %4.1f%% vs %4.1f%%  RR %.3f (%.3f-%.3f)  P=%.3g%s",
          label, df$Outcome, df$ExposedN, df$ComparatorN, 100 * df$ExposedRisk,
          100 * df$ComparatorRisk, df$RR, df$RR_Lower, df$RR_Upper, df$PValue,
          ifelse(is.na(df$NNT), "", sprintf("  NNT %d (rounded up from %.1f)", ceiling(df$NNT), df$NNT)))
}
message(paste(fmt(each_statin, each_statin$Statin), collapse = "\n"))
message(paste(fmt(active, "Atorva vs other"), collapse = "\n"))
message_done("TriNetX outcome tables written to ", PATHS$tables)
