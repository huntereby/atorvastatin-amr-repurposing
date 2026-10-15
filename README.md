# Atorvastatin as a repurposed therapy for kidney transplant antibody-mediated rejection

Code and data for:

> Eby H, Ryan W V, Wood JF, Davuluri T, Vergis J, Arquette M, Ratna A, Rees M, Weistuch C, McCullumsmith R.
> **Insights into renal transplant rejection: potential role of atorvastatin as a repurposed therapy.**
> Submitted to the *American Journal of Transplantation* (2026).

A whole-blood antibody-mediated rejection (AMR) transcriptional signature was queried against the LINCS L1000 chemical-perturbagen library to find compounds whose expression profiles oppose AMR. The top hit, atorvastatin, was then evaluated in propensity-matched kidney transplant cohorts in TriNetX. Statin structural similarity, predicted protein targets, and LINCS gene-knockdown connectivity were then used to prioritize mechanisms specific to atorvastatin.

## Repository layout

```
├── run_all.sh                  run the whole pipeline in order
├── R/utils.R                   shared paths and helpers
├── scripts/                    numbered pipeline steps (see below)
├── data/
│   ├── raw/
│   │   ├── trinetx/            TriNetX aggregate exports (5 analyses)
│   │   ├── target_prediction/  SwissTargetPrediction and SuperPred results (4 statins)
│   │   ├── lincs/              atorvastatin L1000 signature LINCSCP_50693 (iLINCS)
│   │   └── statin_structures.csv   SMILES and molecular formulas for 8 statins
│   └── cache/ilincs/           iLINCS query results reported in the manuscript
├── results/
│   ├── intermediate/           differential expression and AMR gene set
│   ├── tables/                 every table behind a figure or reported number
│   └── figures/                Figures 1–4 and S1 (PDF, 600-dpi PNG and TIFF)
└── env/                        R install script, package versions, Python requirements
```

## Pipeline

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 01 | `01_amr_signature.R` | Downloads GSE50084 (GPL6244) from GEO. Runs limma for whole-blood AMR (DSA+, biopsy-proven; n = 28) vs DSA-negative (n = 12), with BH FDR. Defines the AMR gene set: genes in the top 10% by log2FC in either direction that also have FDR < 0.05. | `dge_blood_amr_vs_dsa_negative.csv`, `amr_gene_set.csv` (2,046 genes) |
| 02 | `02_lincs_drug_screen.R` | drugfindR/iLINCS chemical-perturbagen query: top/bottom 10% of L1000 landmark genes, paired up/down. Ranks discordant compounds, filters to FDA-approved drugs (Drugs@FDA), and annotates mechanisms of action (Broad Repurposing Hub, 2025-08-18). | `lincs_top10_fda_approved.csv`, `lincs_top100_moa_counts.csv` |
| 03 | `03_trinetx_outcomes.R` | Parses the TriNetX outcome exports for rejection (T86.11) and graft failure (T86.12), days 1–365. Reports risk ratios, risk differences and NNT. | `trinetx_each_statin_vs_no_statin.csv`, `trinetx_atorvastatin_vs_other_statins.csv` |
| 04 | `04a_table_s1_four_statins.py`, `04b_table_s2_atorvastatin_vs_other_statins.py` | Baseline characteristics before and after propensity matching. | `Table_S1_*.docx/csv`, `Table_S2_*.docx/csv` |
| 05 | `05_structural_similarity.py` | Checks each SMILES against its molecular formula. Computes RDKit Morgan fingerprints (r = 2, 2,048 bits; ECFP4) and Tanimoto similarity for all 28 statin pairs. | `statin_tanimoto_matrix.csv`, `statin_pairwise_tanimoto.csv` |
| 06 | `06_target_prediction.R` | Takes the union of SwissTargetPrediction (top 100) and SuperPred targets for each statin, maps UniProt IDs to HGNC symbols, assigns targets to the 15 four-set Venn regions, and intersects them with the AMR gene set. | `target_union_membership.csv`, `target_venn_regions_amr.csv` |
| 07 | `07_kd_connectivity.R` | Queries the LINCS knockdown library with the AMR signature and with the atorvastatin signature (top/bottom 100 landmark genes each). Scores the atorvastatin-only and atorvastatin + rosuvastatin-only targets. | `kd_connectivity_*_targets.csv` |
| 08 | `08_main_figures.R` | Figures 1–4 from the tables above | `Figure1`–`Figure4` |
| 09 | `09_figure_s1.R` | Figure S1 (volcano plot and AMR gene-set construction) | `FigureS1` |

## Manuscript item → source

| Manuscript item | Script | Table(s) in `results/tables/` |
|---|---|---|
| Figure 1 | 02, 08 | `lincs_top10_fda_approved.csv` |
| Figure 2A | 03, 08 | `trinetx_each_statin_vs_no_statin.csv` |
| Figure 2B | 03, 08 | `trinetx_atorvastatin_vs_other_statins.csv` |
| Figure 3A–B | 05, 08 | `statin_tanimoto_matrix.csv`, `statin_pairwise_tanimoto.csv`, `statin_family_tanimoto_summary.csv` |
| Figure 4A | 06, 08 | `target_union_membership.csv`, `target_venn_regions_amr.csv` |
| Figure 4B | 07, 08 | `kd_connectivity_atorvastatin_only_targets.csv`, `kd_connectivity_atorvastatin_rosuvastatin_targets.csv` |
| Figure S1 | 01, 09 | `figure_s1_gene_counts.csv` |
| Tables S1A/B | 04a | `Table_S1_four_statins_vs_no_statin.docx` |
| Table S2 | 04b | `Table_S2_atorvastatin_vs_other_statins.docx` |
| Screen summary (7,301 compounds; 4,932 discordant) | 02 | `results/intermediate/lincs_cp_results.csv` |

## Running it

Requirements: R ≥ 4.4 with Bioconductor, and Python ≥ 3.9. Package versions used for the manuscript are in `env/r_package_versions.csv` and `env/requirements.txt`.

```bash
Rscript env/install_r_packages.R
python3 -m pip install -r env/requirements.txt
bash run_all.sh
```

Run every script from the repository root. Each step reads only files in `data/` or outputs of earlier steps, so you can rerun any single step, e.g. `Rscript scripts/08_main_figures.R`. Steps 01 and 02 need internet access (GEO, Drugs@FDA, Broad Repurposing Hub).

### iLINCS queries are cached

The LINCS connectivity results in `data/cache/ilincs/` are the queries reported in the manuscript. By default, steps 02 and 07 read these files instead of calling the iLINCS API, because the LINCS libraries and the iLINCS service change over time. To re-query, run `REFRESH_ILINCS=1 bash run_all.sh`; the new results may differ slightly from the published ones.

| Cache file | Query |
|---|---|
| `cp_query_amr_signature.csv` | AMR signature vs chemical perturbagens (CP), \|similarity\| ≥ 0.2 |
| `kd_query_amr_signature.csv` | AMR signature vs consensus gene knockdowns (KD) |
| `kd_query_atorvastatin_signature.csv` | Atorvastatin LINCSCP_50693 (PL21, 10 µM, 6 h) vs KD |

iLINCS returned only signatures with \|similarity\| ≥ 0.2 for both KD queries, so targets whose knockdown connectivity is below that on either axis have no score.

## Data sources

| Source | Access |
|---|---|
| GSE50084 (Ó Broin et al., *Genom Data* 2014) | Public, GEO; downloaded by step 01 |
| LINCS L1000 via iLINCS (Pilarczyk et al., *Nat Commun* 2022) | Public; queried with [drugfindR](https://github.com/CogDisResLab/drugfindR) 1.0.0 |
| Drugs@FDA | Public; downloaded by step 02 |
| Broad Drug Repurposing Hub, release 2025-08-18 | Public; downloaded by step 02 |
| SwissTargetPrediction, SuperPred 3.0 | Public web servers; exported results in `data/raw/target_prediction/` |
| TriNetX US Collaborative Network | Aggregate, de-identified exports in `data/raw/trinetx/`. Patient-level data cannot be shared; the cohorts can be rebuilt by TriNetX users from the definitions in the manuscript Methods and in each export's `_Read_Me.txt`. |

## TriNetX cohorts

Each folder in `data/raw/trinetx/` is a TriNetX "Compare Outcomes" export, with 1:1 propensity-score matching. Index event: renal allotransplantation (CPT 1008109). Outcome 1 is kidney transplant rejection (T86.11); outcome 2 is kidney transplant failure (T86.12); both are assessed on days 1–365. TriNetX labels both as "Unnamed Outcome" in its metadata.

| Folder | Comparison | Run on | HCOs online |
|---|---|---|---|
| `atorvastatin_vs_no_statin` | Atorvastatin vs statin-naïve | 2026-09-17 | 75 of 75 |
| `rosuvastatin_vs_no_statin` | Rosuvastatin vs statin-naïve | 2026-09-16 | 75 of 75 |
| `simvastatin_vs_no_statin` | Simvastatin vs statin-naïve | 2026-09-17 | 75 of 75 |
| `pravastatin_vs_no_statin` | Pravastatin vs statin-naïve | 2026-09-17 | 75 of 75 |
| `atorvastatin_vs_other_statins` | Atorvastatin vs fluvastatin, rosuvastatin, pitavastatin, pravastatin or simvastatin (active comparator) | 2026-09-30 | 74 of 75 |

In every analysis, patients with the outcome before the time window are kept, not excluded. NNTs in the result tables are unrounded; the manuscript reports them rounded up to the next whole patient.

## License

Code is released under the MIT License (see `LICENSE`). Third-party data keep their original terms.

## Contact

Hunter Eby, Department of Neurosciences, University of Toledo College of Medicine and Life Sciences: hunter.eby@rockets.utoledo.edu
