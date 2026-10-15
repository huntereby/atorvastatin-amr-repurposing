# Install the R packages used by the pipeline.
# Versions used for the manuscript are listed in env/r_package_versions.csv
# (R 4.6.1, Bioconductor 3.23).

cran <- c("dplyr", "readr", "tidyr", "stringr", "purrr", "tibble", "jsonlite",
          "ggplot2", "ggrepel", "ggVennDiagram", "patchwork", "scales")
bioc <- c("GEOquery", "Biobase", "limma", "AnnotationDbi", "org.Hs.eg.db")

repos <- "https://cloud.r-project.org"
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = repos)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

# drugfindR (iLINCS client; CogDisResLab). Try Bioconductor first, then GitHub.
if (!requireNamespace("drugfindR", quietly = TRUE)) {
  ok <- tryCatch({ BiocManager::install("drugfindR", update = FALSE, ask = FALSE); TRUE },
                 error = function(e) FALSE)
  if (!ok || !requireNamespace("drugfindR", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
    remotes::install_github("CogDisResLab/drugfindR")
  }
}
