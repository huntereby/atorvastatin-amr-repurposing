#!/usr/bin/env bash
# Run the full pipeline from the repository root.
#
#   bash run_all.sh            # use cached iLINCS queries (reproduces the manuscript)
#   REFRESH_ILINCS=1 bash run_all.sh   # re-query iLINCS (results may differ)
#
# Step 01 downloads GSE50084 from GEO; step 02 downloads Drugs@FDA and the
# Broad Drug Repurposing Hub annotation into data/cache/external/.
set -euo pipefail
cd "$(dirname "$0")"
PYTHON="${PYTHON:-python3}"

Rscript scripts/01_amr_signature.R
Rscript scripts/02_lincs_drug_screen.R
Rscript scripts/03_trinetx_outcomes.R
"$PYTHON" scripts/04a_table_s1_four_statins.py
"$PYTHON" scripts/04b_table_s2_atorvastatin_vs_other_statins.py
"$PYTHON" scripts/05_structural_similarity.py
Rscript scripts/06_target_prediction.R
Rscript scripts/07_kd_connectivity.R
Rscript scripts/08_main_figures.R
Rscript scripts/09_figure_s1.R
echo "Pipeline complete. Tables: results/tables/  Figures: results/figures/"
