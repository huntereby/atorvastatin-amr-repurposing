#!/usr/bin/env Rscript
# Step 09 - Figure S1: whole-blood differential expression (GSE50084) and
# construction of the AMR gene set used for target overlap.
#
# Panel A: volcano plot. Genes at FDR < 0.05 are coloured by direction; those
#          also in the top 10% by fold change in each direction (the AMR gene
#          set) are shown darker. Dotted lines mark the top-10% fold-change
#          cutoffs.
# Panel B: gene counts at each filtering step, split by direction.
#
# Inputs:  results/intermediate/dge_blood_amr_vs_dsa_negative.csv, amr_top10pct_{up,down}.csv
# Outputs: results/figures/FigureS1.{pdf,png,tiff}
#          results/tables/figure_s1_gene_counts.csv

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(tidyr)
})

dge_file <- file.path(PATHS$intermediate, "dge_blood_amr_vs_dsa_negative.csv")
up_file <- file.path(PATHS$intermediate, "amr_top10pct_up.csv")
down_file <- file.path(PATHS$intermediate, "amr_top10pct_down.csv")
stopifnot(file.exists(dge_file), file.exists(up_file), file.exists(down_file))

fdr_cutoff <- 0.05
n_labels <- 8L

dge <- read_csv(dge_file, show_col_types = FALSE)
top_up <- read_csv(up_file, show_col_types = FALSE)
top_down <- read_csv(down_file, show_col_types = FALSE)

up_fc_cutoff <- min(top_up$logFC)
down_fc_cutoff <- max(top_down$logFC)

class_levels <- c(
  "Not significant",
  "FDR < 0.05, higher in AMR",
  "FDR < 0.05, lower in AMR",
  "Target-overlap set, higher in AMR",
  "Target-overlap set, lower in AMR"
)

dge <- dge |>
  mutate(
    Significant = adj.P.Val < fdr_cutoff,
    InTop10 = GeneSymbol %in% c(top_up$GeneSymbol, top_down$GeneSymbol),
    Class = case_when(
      !Significant ~ class_levels[1],
      InTop10 & logFC > 0 ~ class_levels[4],
      InTop10 & logFC < 0 ~ class_levels[5],
      logFC > 0 ~ class_levels[2],
      TRUE ~ class_levels[3]
    ),
    Class = factor(Class, levels = class_levels),
    NegLog10FDR = -log10(adj.P.Val)
  )

counts <- list(
  tested = nrow(dge),
  sig_up = sum(dge$Significant & dge$logFC > 0),
  sig_down = sum(dge$Significant & dge$logFC < 0),
  set_up = sum(dge$Class == class_levels[4]),
  set_down = sum(dge$Class == class_levels[5]),
  top_n = nrow(top_up)
)
stopifnot(nrow(top_up) == nrow(top_down))

class_colours <- c(
  "Not significant" = "grey82",
  "FDR < 0.05, higher in AMR" = "#F4A582",
  "FDR < 0.05, lower in AMR" = "#92C5DE",
  "Target-overlap set, higher in AMR" = "#B2182B",
  "Target-overlap set, lower in AMR" = "#2166AC"
)

theme_supp <- function() {
  theme_classic(base_size = 10, base_family = "Helvetica") +
    theme(
      plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
      plot.subtitle = element_text(size = 8.8, colour = "grey25", hjust = 0),
      plot.tag = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.margin = margin(8, 14, 8, 8)
    )
}

# ---------------------------------------------------------------------------
# Panel A: volcano
# ---------------------------------------------------------------------------
label_genes <- dge |>
  filter(
    Class %in% class_levels[4:5],
    !grepl("^(RNA5SP|RN7S|RNU|SNOR|MIR|LINC|LOC)|-AS[0-9]*$|P[0-9]+$", GeneSymbol)
  ) |>
  group_by(Direction = sign(logFC)) |>
  slice_max(NegLog10FDR * abs(logFC), n = n_labels / 2, with_ties = FALSE) |>
  ungroup()

x_lim <- max(abs(dge$logFC)) * 1.05

panel_a <- ggplot(dge, aes(logFC, NegLog10FDR, colour = Class)) +
  geom_point(data = filter(dge, Class == class_levels[1]), size = 0.5, alpha = 0.5) +
  geom_point(data = filter(dge, Class %in% class_levels[2:3]), size = 0.8, alpha = 0.8) +
  geom_point(data = filter(dge, Class %in% class_levels[4:5]), size = 0.9, alpha = 0.85) +
  geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed",
             colour = "grey30", linewidth = 0.35) +
  geom_vline(xintercept = c(down_fc_cutoff, up_fc_cutoff), linetype = "dotted",
             colour = "grey30", linewidth = 0.4) +
  annotate("text", x = x_lim, y = -log10(fdr_cutoff), label = "FDR = 0.05",
           hjust = 1, vjust = -0.5, size = 2.7, colour = "grey30", family = "Helvetica") +
  geom_text_repel(
    data = label_genes, aes(label = GeneSymbol),
    size = 2.7, fontface = "italic", family = "Helvetica", colour = "black",
    min.segment.length = 0, segment.colour = "grey45", segment.size = 0.3,
    box.padding = 0.35, max.overlaps = Inf, seed = 1, show.legend = FALSE
  ) +
  scale_colour_manual(values = class_colours, name = NULL, drop = FALSE) +
  scale_x_continuous(limits = c(-x_lim, x_lim)) +
  labs(
    title = "Differential expression, whole blood",
    subtitle = sprintf(
      "Dotted lines: top-10%% fold-change cutoffs (log2FC %.2f and +%.2f)",
      down_fc_cutoff, up_fc_cutoff
    ),
    x = expression(bold(log[2] * " fold change (AMR vs DSA-negative)")),
    y = expression(bold(-log[10] * " FDR"))
  ) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2.6, alpha = 1))) +
  theme_supp() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

# ---------------------------------------------------------------------------
# Panel B: gene counts at each filtering step
# ---------------------------------------------------------------------------
step_levels <- c(
  "Genes tested",
  "FDR < 0.05",
  "FDR < 0.05 and\ntop 10% fold change\n(target-overlap set)"
)

funnel <- tibble(
  Step = factor(rep(step_levels, each = 2), levels = rev(step_levels)),
  Direction = rep(c("Higher in AMR", "Lower in AMR"), times = 3),
  Genes = c(NA, NA, counts$sig_up, counts$sig_down, counts$set_up, counts$set_down)
)
totals <- tibble(
  Step = factor(step_levels, levels = rev(step_levels)),
  Total = c(counts$tested, counts$sig_up + counts$sig_down, counts$set_up + counts$set_down)
)

direction_colours <- c("Higher in AMR" = "#B2182B", "Lower in AMR" = "#2166AC")

panel_b <- ggplot(filter(funnel, !is.na(Genes)), aes(Genes, Step, fill = Direction)) +
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(
    aes(label = format(Genes, big.mark = ",", trim = TRUE)),
    position = position_stack(vjust = 0.5, reverse = TRUE),
    colour = "white", size = 3.1, fontface = "bold", family = "Helvetica"
  ) +
  geom_text(
    data = totals, aes(x = pmin(Total, 2900), y = Step, label = paste0("n = ", format(Total, big.mark = ",", trim = TRUE))),
    inherit.aes = FALSE, hjust = -0.12, size = 3.2, fontface = "bold", family = "Helvetica"
  ) +
  annotate(
    "segment", x = 0, xend = 2900, y = 3, yend = 3,
    colour = "grey55", linewidth = 6, lineend = "butt"
  ) +
  annotate(
    "text", x = 1450, y = 3, label = "all annotated genes (axis truncated)",
    colour = "white", size = 2.8, fontface = "bold", family = "Helvetica"
  ) +
  scale_fill_manual(values = direction_colours, name = NULL) +
  scale_x_continuous(
    limits = c(0, 3600), breaks = seq(0, 2500, by = 500),
    labels = scales::label_comma(), expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "AMR gene set construction",
    subtitle = sprintf(
      "Top 10%% = %s genes per direction ranked by fold change",
      format(counts$top_n, big.mark = ",")
    ),
    x = "Genes",
    y = NULL
  ) +
  theme_supp() +
  theme(
    legend.position = "bottom",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 9, face = "bold"),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

figure <- (panel_a | panel_b) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Whole-blood transcriptome in AMR versus DSA-negative recipients (GSE50084)",
    subtitle = "DSA rejection (n = 28) vs No DSA (n = 12); limma, Benjamini-Hochberg FDR",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, family = "Helvetica"),
      plot.subtitle = element_text(size = 9.5, colour = "grey25", family = "Helvetica")
    )
  )

write_csv(
  tibble(
    Metric = c("Genes tested", "FDR < 0.05 higher in AMR", "FDR < 0.05 lower in AMR",
               "Top-10% genes per direction", "Up log2FC cutoff", "Down log2FC cutoff",
               "Target-overlap set higher in AMR", "Target-overlap set lower in AMR"),
    Value = c(counts$tested, counts$sig_up, counts$sig_down, counts$top_n,
              round(up_fc_cutoff, 4), round(down_fc_cutoff, 4),
              counts$set_up, counts$set_down)
  ),
  file.path(PATHS$tables, "figure_s1_gene_counts.csv")
)

save_figure <- function(plot, base, width, height) {
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 400, bg = "white")
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height,
                 family = "Helvetica", encoding = "WinAnsi.enc", useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()

  # The macOS Quartz TIFF device does not implement compression. Render a
  # lossless 600-dpi PNG first, then convert it to an LZW-compressed TIFF.
  high_res_png <- tempfile(pattern = "figure_s1_", fileext = ".png")
  ggsave(high_res_png, plot, width = width, height = height, dpi = 600, bg = "white")
  if (!file.exists("/usr/bin/sips")) {
    ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = 600,
           bg = "white", compression = "lzw")
    return(invisible())
  }
  system2(
    "/usr/bin/sips",
    c("-s", "format", "tiff", "-s", "formatOptions", "lzw",
      "-s", "dpiWidth", "600", "-s", "dpiHeight", "600",
      high_res_png, "--out", paste0(base, ".tiff")),
    stdout = FALSE, stderr = FALSE
  )
  unlink(high_res_png)
}

save_figure(figure, file.path(PATHS$figures, "FigureS1"), 13.5, 7.2)

message(sprintf(
  "Tested %d; FDR<0.05 %d up / %d down; target-overlap set %d up / %d down",
  counts$tested, counts$sig_up, counts$sig_down, counts$set_up, counts$set_down
))
