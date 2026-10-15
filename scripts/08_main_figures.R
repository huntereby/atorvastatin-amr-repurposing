#!/usr/bin/env Rscript
# Step 08 - Main manuscript figures (Figures 1-4), drawn from the tables
# written by steps 02-07. Nothing upstream is recomputed.
#
# One shared print style: 6.5 in wide (3.5 in for Figure 1), Helvetica 8 pt,
# and one colour scheme used with the same meaning everywhere:
#   blue = favourable direction (lower risk; target opposes the AMR signature)
#   red  = unfavourable direction (higher risk)
#   grey = not significant / neutral
#
# Outputs: results/figures/Figure{1,2,3,4}.{pdf,png,tiff}

source(file.path("R", "utils.R"))
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggVennDiagram)
  library(patchwork)
  library(tidyr)
})

# ---------------------------------------------------------------------------
# Shared style
# ---------------------------------------------------------------------------
FONT <- "Helvetica"
FONT_PT <- 8
TEXT_SIZE <- FONT_PT / .pt # geom_text size giving 8 pt
FIG_WIDTH <- 6.5

BLUE <- "#2166AC"
BLUE_LIGHT <- "#92C5DE"
RED <- "#B2182B"
GREY <- "#8C8C8C"
GREY_LIGHT <- "#D9D9D9"
NAVY <- "#08306B"
INK <- "#000000"
INK_MUTED <- "#4D4D4D"

theme_print <- function() {
  theme_classic(base_size = FONT_PT, base_family = FONT) +
    theme(
      text = element_text(size = FONT_PT, colour = INK),
      plot.title = element_text(size = FONT_PT, face = "bold", hjust = 0, margin = margin(b = 4)),
      plot.subtitle = element_text(size = FONT_PT, colour = INK_MUTED, hjust = 0),
      plot.tag = element_text(size = FONT_PT, face = "bold"),
      axis.text = element_text(size = FONT_PT, colour = INK),
      axis.title = element_text(size = FONT_PT, colour = INK),
      axis.line = element_line(colour = INK, linewidth = 0.3),
      axis.ticks = element_line(colour = INK, linewidth = 0.3),
      legend.text = element_text(size = FONT_PT),
      legend.title = element_text(size = FONT_PT),
      legend.key.size = unit(3, "mm"),
      strip.text = element_text(size = FONT_PT, face = "bold"),
      strip.background = element_rect(fill = "#F0F0F0", colour = NA),
      plot.margin = margin(4, 4, 4, 4)
    )
}

annotate_tags <- function(fig) {
  fig + plot_annotation(tag_levels = "A", theme = theme(text = element_text(family = FONT)))
}

# PDF (vector), 600-dpi PNG, and 600-dpi LZW-compressed TIFF for submission.
save_print <- function(plot, name, height, width = FIG_WIDTH) {
  base <- file.path(PATHS$figures, name)
  grDevices::pdf(paste0(base, ".pdf"), width = width, height = height, family = FONT,
                 encoding = "WinAnsi.enc", useDingbats = FALSE) # WinAnsi keeps the en dash
  print(plot)
  grDevices::dev.off()
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 600, bg = "white")

  if (file.exists("/usr/bin/sips")) {
    # The macOS Quartz TIFF device does not compress; convert the PNG instead.
    system2("/usr/bin/sips",
            c("-s", "format", "tiff", "-s", "formatOptions", "lzw",
              "-s", "dpiWidth", "600", "-s", "dpiHeight", "600",
              paste0(base, ".png"), "--out", paste0(base, ".tiff")),
            stdout = FALSE, stderr = FALSE)
  } else {
    ggsave(paste0(base, ".tiff"), plot, width = width, height = height, dpi = 600,
           bg = "white", compression = "lzw")
  }
  message("Wrote ", base, ".{pdf,png,tiff}")
}

# ---------------------------------------------------------------------------
# Figure 1: FDA-approved drugs most discordant with the AMR signature
# ---------------------------------------------------------------------------
top_drugs <- read_csv(file.path(PATHS$tables, "lincs_top10_fda_approved.csv"), show_col_types = FALSE) |>
  mutate(
    Drug = tools::toTitleCase(tolower(Target)),
    Drug = case_when(Drug == "Amidate" ~ "Etomidate (Amidate)",
                     Drug == "Quinidine Gluconate" ~ "Quinidine gluconate",
                     TRUE ~ Drug),
    Label = sprintf("%s [%d]", Drug, OverallDiscordantRank),
    Label = factor(Label, levels = rev(Label)),
    Highlight = if_else(tolower(Target) == "atorvastatin", "Atorvastatin", "Other")
  )

fig1 <- ggplot(top_drugs, aes(Similarity, Label, fill = Highlight)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Similarity), colour = Highlight), hjust = -0.15,
            size = TEXT_SIZE, family = FONT) +
  scale_colour_manual(values = c(Atorvastatin = "white", Other = INK), guide = "none") +
  geom_vline(xintercept = -0.2, linetype = "dashed", colour = INK_MUTED, linewidth = 0.35) +
  scale_fill_manual(values = c(Atorvastatin = BLUE, Other = BLUE_LIGHT), guide = "none") +
  scale_x_continuous(limits = c(-1, 0), breaks = seq(-1, 0, 0.25),
                     labels = c("-1", "-0.75", "-0.5", "-0.25", "0"), expand = c(0, 0)) +
  labs(title = "Top FDA-approved discordant drugs",
       x = "Connectivity with AMR signature (similarity)", y = NULL) +
  theme_print() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.x = element_line(colour = "#EBEBEB", linewidth = 0.3),
        plot.title.position = "plot")

save_print(fig1, "Figure1", 3.0, width = 3.5)

# ---------------------------------------------------------------------------
# Figure 2: TriNetX forest plots
# ---------------------------------------------------------------------------
outcome_levels <- unname(OUTCOME_LABELS)
direction_colours <- c("Lower risk" = BLUE, "Higher risk" = RED, "Not significant" = GREY)
rr_breaks <- c(0.5, 0.75, 1, 1.5, 2)
label_x <- 2.35 # estimate column, just right of the plotting area

forest_panel <- function(data, y_var, title) {
  ggplot(data, aes(RR, .data[[y_var]], colour = Direction)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = INK_MUTED, linewidth = 0.35) +
    geom_errorbarh(aes(xmin = RR_Lower, xmax = RR_Upper), height = 0.25, linewidth = 0.6) +
    geom_point(size = 2.2) +
    geom_text(aes(x = label_x, label = EstimateLabel), hjust = 0, colour = INK,
              size = TEXT_SIZE, family = FONT) +
    scale_colour_manual(values = direction_colours, limits = names(direction_colours), name = NULL, drop = FALSE) +
    scale_x_log10(breaks = rr_breaks, labels = c("0.5", "0.75", "1", "1.5", "2")) +
    coord_cartesian(xlim = c(0.5, 2), clip = "off") +
    labs(title = title, x = "Risk ratio (95% CI), log scale", y = NULL) +
    theme_print() +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_line(colour = "#EBEBEB", linewidth = 0.3),
      plot.margin = margin(4, 92, 4, 4)
    )
}

prep_forest <- function(df) {
  df |>
    mutate(
      Direction = factor(Direction, levels = names(direction_colours)),
      EstimateLabel = sprintf("%.2f (%.2f–%.2f)", RR, RR_Lower, RR_Upper)
    )
}

forest_each <- read_csv(file.path(PATHS$tables, "trinetx_each_statin_vs_no_statin.csv"), show_col_types = FALSE) |>
  prep_forest() |>
  mutate(
    Outcome = factor(Outcome, levels = outcome_levels),
    Statin = factor(Statin, levels = rev(STATINS_CLINICAL))
  )
forest_active <- read_csv(file.path(PATHS$tables, "trinetx_atorvastatin_vs_other_statins.csv"), show_col_types = FALSE) |>
  prep_forest() |>
  mutate(Outcome = factor(Outcome, levels = rev(outcome_levels)))

fig2_a <- forest_panel(forest_each, "Statin", "Each statin vs matched no-statin controls") +
  facet_wrap(~Outcome, ncol = 1, scales = "free_y") +
  theme(axis.title.x = element_blank())
fig2_b <- forest_panel(forest_active, "Outcome", "Atorvastatin vs other statins") +
  guides(colour = "none") # one shared key, drawn from panel A

fig2 <- annotate_tags(
  (fig2_a / fig2_b) +
    plot_layout(heights = c(8.6, 2), guides = "collect") &
    theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
)
save_print(fig2, "Figure2", 5.2)

# ---------------------------------------------------------------------------
# Figure 3: statin structural similarity
# ---------------------------------------------------------------------------
tanimoto <- read.csv(file.path(PATHS$tables, "statin_tanimoto_matrix.csv"), row.names = 1, check.names = FALSE)
statin_order <- colnames(tanimoto)
type_ii <- c("Atorvastatin", "Rosuvastatin", "Fluvastatin", "Pitavastatin", "Cerivastatin")
type_i <- c("Simvastatin", "Lovastatin", "Pravastatin")
stopifnot(identical(statin_order, c(type_ii, type_i)))

matrix_long <- as.data.frame(as.table(as.matrix(tanimoto))) |>
  setNames(c("drug_y", "drug_x", "tanimoto")) |>
  mutate(drug_x = factor(drug_x, levels = statin_order),
         drug_y = factor(drug_y, levels = rev(statin_order)))
n_ii <- length(type_ii)
n_all <- length(statin_order)

fig3_a <- ggplot(matrix_long, aes(drug_x, drug_y, fill = tanimoto)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", tanimoto), colour = tanimoto > 0.5),
            size = TEXT_SIZE, family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK)) +
  scale_fill_gradientn(colours = c("#FFFFFF", "#C6DBEF", BLUE, NAVY), limits = c(0, 1),
                       breaks = c(0, 0.5, 1), name = "Tanimoto") +
  geom_vline(xintercept = n_ii + 0.5, colour = INK, linewidth = 0.6) +
  geom_hline(yintercept = length(type_i) + 0.5, colour = INK, linewidth = 0.6) +
  annotate("text", x = (1 + n_ii) / 2, y = n_all + 0.9, label = "Type II (synthetic)",
           size = TEXT_SIZE, family = FONT, fontface = "bold") +
  annotate("text", x = (n_ii + 1 + n_all) / 2, y = n_all + 0.9, label = "Type I (fungal)",
           size = TEXT_SIZE, family = FONT, fontface = "bold") +
  coord_equal(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_print() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.key.height = unit(6, "mm"),
    legend.key.width = unit(2.5, "mm"),
    plot.margin = margin(14, 4, 4, 4)
  )

pairwise <- read_csv(file.path(PATHS$tables, "statin_pairwise_tanimoto.csv"), show_col_types = FALSE)
group_levels <- c("Within Type I", "Within Type II", "Between families")
group_n <- table(factor(pairwise$comparison_group, levels = group_levels))
group_labels <- setNames(
  sprintf("%s\n%s\n(n = %d)", c("Within", "Within", "Between"),
          c("type I", "type II", "families"), as.integer(group_n)),
  group_levels
)
pairwise <- pairwise |> mutate(comparison_group = factor(comparison_group, levels = group_levels))

fig3_b <- ggplot(pairwise, aes(comparison_group, tanimoto, fill = comparison_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.35, colour = INK) +
  geom_point(position = position_jitter(width = 0.08, seed = 1), size = 1.4,
             shape = 21, fill = "white", colour = INK, stroke = 0.4) +
  scale_fill_manual(values = c("Within Type I" = BLUE_LIGHT, "Within Type II" = BLUE_LIGHT,
                               "Between families" = GREY_LIGHT), guide = "none") +
  scale_x_discrete(labels = group_labels) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Tanimoto coefficient") +
  theme_print()

save_print(annotate_tags((fig3_a | fig3_b) + plot_layout(widths = c(1.45, 1))), "Figure3", 3.3)

# ---------------------------------------------------------------------------
# Figure 4: predicted targets (Venn) + knockdown connectivity
# ---------------------------------------------------------------------------
membership <- read_csv(file.path(PATHS$tables, "target_union_membership.csv"), show_col_types = FALSE)
gene_lists <- setNames(lapply(STATINS_CLINICAL, function(d) membership$TargetGene[membership[[d]] == 1]),
                       STATINS_CLINICAL)

venn_data <- process_data(Venn(gene_lists))
region_shape <- ggVennDiagram:::get_shape_regionedge(venn_data)
region_label_pos <- ggVennDiagram:::get_shape_regionlabel(venn_data) # includes region counts
set_edge <- ggVennDiagram:::get_shape_setedge(venn_data, color = INK, linetype = "solid", linewidth = 0.5)
set_label <- ggVennDiagram:::get_shape_setlabel(venn_data, size = TEXT_SIZE, color = INK)

fig4_a <- ggplot(mapping = aes(X, Y)) +
  geom_polygon(data = region_shape, aes(group = id), fill = "white", colour = GREY_LIGHT, linewidth = 0.25) +
  geom_path(data = set_edge, aes(group = id), colour = INK, linewidth = 0.5) +
  geom_text(data = set_label, aes(label = name), size = TEXT_SIZE, family = FONT, fontface = "bold") +
  geom_text(data = region_label_pos, aes(label = count), size = TEXT_SIZE, family = FONT) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_y_continuous(expand = expansion(mult = 0.02)) +
  coord_equal(clip = "off") +
  theme_void(base_size = FONT_PT, base_family = FONT) +
  theme(plot.tag = element_text(size = FONT_PT, face = "bold"), plot.margin = margin(4, 4, 4, 4))

# Atorvastatin-only and atorvastatin + rosuvastatin-only targets with KD
# signatures; AMR gene-set members are highlighted and labelled.
target_set_levels <- c("Atorvastatin only", "Atorvastatin + rosuvastatin")
kd <- bind_rows(
  read_csv(file.path(PATHS$tables, "kd_connectivity_atorvastatin_only_targets.csv"), show_col_types = FALSE) |>
    mutate(TargetSet = target_set_levels[1]),
  read_csv(file.path(PATHS$tables, "kd_connectivity_atorvastatin_rosuvastatin_targets.csv"), show_col_types = FALSE) |>
    mutate(TargetSet = target_set_levels[2])
) |>
  mutate(
    AMR = factor(if_else(AMRGene, "AMR signature gene", "Other predicted target"),
                 levels = c("AMR signature gene", "Other predicted target")),
    TargetSet = factor(TargetSet, levels = target_set_levels)
  )
set_counts <- table(kd$TargetSet)
levels(kd$TargetSet) <- sprintf("%s (n = %d)", target_set_levels, as.integer(set_counts))

fig4_b <- ggplot(kd, aes(r_disease, r_ator)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = BLUE, alpha = 0.06) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = BLUE, alpha = 0.06) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = INK_MUTED) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = INK_MUTED) +
  geom_point(data = filter(kd, AMR == "Other predicted target"),
             aes(shape = TargetSet, fill = AMR), size = 1.8, colour = "white", stroke = 0.3) +
  geom_point(data = filter(kd, AMR == "AMR signature gene"),
             aes(shape = TargetSet, fill = AMR), size = 2.6, colour = "white", stroke = 0.4) +
  geom_text_repel(data = filter(kd, AMR == "AMR signature gene"), aes(label = Target),
                  size = TEXT_SIZE, family = FONT, fontface = "bold", colour = INK,
                  segment.size = 0.3, segment.colour = INK_MUTED, min.segment.length = 0,
                  box.padding = 0.6, point.padding = 0.3, seed = 3) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_fill_manual(values = c("AMR signature gene" = BLUE, "Other predicted target" = GREY_LIGHT),
                    breaks = c("AMR signature gene", "Other predicted target"), name = NULL) +
  scale_x_continuous(limits = c(-0.5, 0.5), breaks = seq(-0.4, 0.4, 0.2)) +
  scale_y_continuous(limits = c(-0.55, 0.6), breaks = seq(-0.4, 0.6, 0.2)) +
  guides(
    fill = guide_legend(order = 1, override.aes = list(shape = 21, size = 2.4, colour = "white")),
    shape = guide_legend(order = 2, override.aes = list(fill = INK_MUTED, size = 2.2, colour = "white"))
  ) +
  labs(
    x = expression("KD connectivity with AMR signature (" * italic(r)[disease] * ")"),
    y = expression("KD connectivity with atorvastatin (" * italic(r)[atorvastatin] * ")")
  ) +
  theme_print() +
  theme(legend.position = "bottom")

fig4 <- annotate_tags(
  (fig4_a | fig4_b) + plot_layout(widths = c(1, 1), guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal",
          legend.margin = margin(0, 6, 0, 6), legend.text = element_text(size = FONT_PT))
)
save_print(fig4, "Figure4", 3.9)

message_done("Figures 1-4 written to ", PATHS$figures)
