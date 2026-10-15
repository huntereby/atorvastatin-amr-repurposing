# Shared paths and helpers for the analysis scripts.
# Every script is run from the repository root, e.g. `Rscript scripts/01_amr_signature.R`.

suppressPackageStartupMessages({
  library(readr)
})

PATHS <- list(
  raw_trinetx      = file.path("data", "raw", "trinetx"),
  raw_targets      = file.path("data", "raw", "target_prediction"),
  raw_lincs        = file.path("data", "raw", "lincs"),
  ilincs_cache     = file.path("data", "cache", "ilincs"),
  external_cache   = file.path("data", "cache", "external"),
  intermediate     = file.path("results", "intermediate"),
  tables           = file.path("results", "tables"),
  figures          = file.path("results", "figures")
)

for (p in PATHS) dir.create(p, recursive = TRUE, showWarnings = FALSE)

STATINS_CLINICAL <- c("Atorvastatin", "Rosuvastatin", "Simvastatin", "Pravastatin")

OUTCOME_LABELS <- c(
  "1" = "Graft rejection / AMR episode", # ICD-10-CM T86.11
  "2" = "Graft failure / loss"           # ICD-10-CM T86.12
)

# Set REFRESH_ILINCS=1 to re-query iLINCS instead of reading data/cache/ilincs.
# iLINCS results can change as the LINCS libraries are updated, so the cached
# queries are the ones reported in the manuscript.
refresh_ilincs <- function() identical(Sys.getenv("REFRESH_ILINCS"), "1")

# ---------------------------------------------------------------------------
# TriNetX exports: a free-text preamble followed by small CSV blocks, each
# introduced by its own header line.
# ---------------------------------------------------------------------------
read_trinetx_block <- function(lines, header_pattern, n_rows) {
  start <- grep(header_pattern, lines)[1]
  if (is.na(start)) stop("TriNetX block not found: ", header_pattern, call. = FALSE)
  block <- lines[start:(start + n_rows)]
  read.csv(text = paste(block, collapse = "\n"), check.names = FALSE, stringsAsFactors = FALSE)
}

read_trinetx_outcome <- function(study_dir, outcome_num) {
  path <- file.path(study_dir, sprintf("Outcome_%d_Result_a_MOA_table.csv", outcome_num))
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  cohorts <- read_trinetx_block(
    lines, "^Cohort,Cohort Name,Patients in Cohort,Patients with Outcome,Risk$", 2
  )
  rd <- read_trinetx_block(lines, "^Risk Difference,95 % CI Lower,95 % CI Upper,z,p$", 1)
  rr <- read_trinetx_block(lines, "^Risk Ratio,95 % CI Lower,95 % CI Upper$", 1)

  data.frame(
    Outcome = OUTCOME_LABELS[[as.character(outcome_num)]],
    ExposedN = cohorts[1, "Patients in Cohort"],
    ExposedEvents = cohorts[1, "Patients with Outcome"],
    ExposedRisk = cohorts[1, "Risk"],
    ComparatorN = cohorts[2, "Patients in Cohort"],
    ComparatorEvents = cohorts[2, "Patients with Outcome"],
    ComparatorRisk = cohorts[2, "Risk"],
    RiskDifference = rd[1, 1],
    RD_Lower = rd[1, 2],
    RD_Upper = rd[1, 3],
    PValue = rd[1, "p"],
    RR = rr[1, 1],
    RR_Lower = rr[1, 2],
    RR_Upper = rr[1, 3],
    stringsAsFactors = FALSE
  )
}

# Number needed to treat = 1 / absolute risk reduction; its CI is the
# reciprocal of the risk-difference CI limits. Only defined when the whole
# risk-difference CI lies below zero (benefit).
add_nnt <- function(df) {
  benefit <- df$RD_Upper < 0
  df$NNT <- ifelse(benefit, 1 / -df$RiskDifference, NA_real_)
  df$NNT_Lower <- ifelse(benefit, 1 / -df$RD_Lower, NA_real_)
  df$NNT_Upper <- ifelse(benefit, 1 / -df$RD_Upper, NA_real_)
  df
}

# ---------------------------------------------------------------------------
# iLINCS connectivity results: one row per LINCS signature. Collapse to one
# value per gene/compound using the median when several signatures exist.
# ---------------------------------------------------------------------------
collapse_kd_by_gene <- function(df) {
  df <- df[is.finite(df$Similarity), ]
  df$Target <- toupper(trimws(df$Target))
  out <- do.call(rbind, lapply(split(df, df$Target), function(x) {
    data.frame(
      Target = x$Target[1],
      Similarity = stats::median(x$Similarity),
      pValue = stats::median(x$pValue),
      n_signatures = nrow(x)
    )
  }))
  rownames(out) <- NULL
  out
}

message_done <- function(...) message("[done] ", ...)
