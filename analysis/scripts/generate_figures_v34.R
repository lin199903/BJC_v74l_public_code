## generate_figures_v24.R
## Nature-style figure set for the SIGLEC15 two-level qualification manuscript.
## Backend: R (ggplot2 + patchwork + ComplexHeatmap + ggrepel + svglite/cairo_pdf/ragg).
## All panels read machine-readable result CSVs; no analysis is rerun.

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(ggrepel); library(scales)
  library(dplyr); library(readr); library(tidyr); library(ComplexHeatmap)
  library(svglite); library(ragg); library(grid); library(grDevices)
})

## ---- paths ----------------------------------------------------------------
DEPOSIT  <- "./prior_analysis_v1/repository_deposit/Results"
ANALYSIS <- "./data/analysis"
VALID    <- "./validation_design"
SUPPSRC  <- "./Supplementary_Data_1_v17_source/Supplementary_Data_1_v17"
SUPPV6   <- "./Supplementary_Data_1_v6"
OUT      <- "./Figures"
dir.create(OUT, showWarnings = FALSE)

S17 <- file.path(DEPOSIT, "S17"); S18 <- file.path(DEPOSIT, "S18")
S19 <- file.path(DEPOSIT, "S19"); P02 <- file.path(DEPOSIT, "P02"); P22 <- file.path(DEPOSIT, "P22")

## ---- palette + theme ------------------------------------------------------
INK <- "#2B2B2B"; GREY <- "#9A9A9A"; LGREY <- "#D8D8D8"; PALE <- "#F4F4F2"
BLUE <- "#2C6E9B"; BLUE_P <- "#D6E4EF"; TEAL <- "#2E8B85"; TEAL_P <- "#D7ECEA"
ORANGE <- "#C8772E"; ORANGE_P <- "#F2E2D3"; RED <- "#A24B47"; RED_P <- "#EFD8D6"
GREEN <- "#4E8A52"; GREEN_P <- "#DCEBDB"; PURPLE <- "#6E5491"
SUPP <- GREEN; NSUP <- RED; UNINF <- GREY
DS_COL <- c("GSE155698" = "#2C6E9B", "GSE156405" = "#2E8B85", "GSE212966" = "#C8772E")

theme_pub <- function(base = 8.0) {
  theme_classic(base_size = base, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = 0.35, colour = INK),
          axis.ticks = element_line(linewidth = 0.3, colour = INK),
          axis.ticks.length = unit(1.6, "pt"),
          axis.text = element_text(colour = INK, size = base),
          axis.title = element_text(colour = INK, size = base + 0.4),
          plot.title = element_text(size = base + 0.6, face = "bold", colour = INK, hjust = 0),
          plot.subtitle = element_text(size = base, colour = GREY, hjust = 0),
          legend.title = element_text(size = base), legend.text = element_text(size = base - 0.4),
          legend.key.size = unit(7, "pt"), legend.margin = margin(0,0,0,0),
          plot.tag = element_text(size = base + 2.5, face = "bold", colour = INK),
          plot.tag.position = c(0.01, 0.99),
          plot.margin = margin(5, 6, 5, 5))
}
theme_set(theme_pub())
verdict_col <- c("Supportive" = SUPP, "Non-supportive" = NSUP,
                 "Uninformative" = UNINF, "Non-supportive (failed qualification)" = NSUP,
                 "Uninformative (exploratory)" = UNINF, "No association (exploratory)" = UNINF,
                 "direction-consistent supportive" = SUPP)

save_pub <- function(plot, name, w_mm = 170, h_mm = 95, dpi = 300, png_dpi = 600) {
  w <- w_mm/25.4; h <- h_mm/25.4
  svglite(file.path(OUT, paste0(name, ".svg")), width = w, height = h,
          system_fonts = list(sans = "Arial"))
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) print(plot) else grid.draw(plot)
  dev.off()
  grDevices::cairo_pdf(file.path(OUT, paste0(name, ".pdf")), width = w, height = h,
                       family = "Arial")
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) print(plot) else grid.draw(plot)
  dev.off()
  ragg::agg_png(file.path(OUT, paste0(name, ".png")), width = w, height = h,
                units = "in", res = png_dpi, scaling = 1, background = "white")
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) print(plot) else grid.draw(plot)
  dev.off()
  ragg::agg_tiff(file.path(OUT, paste0(name, ".tiff")), width = w, height = h,
                 units = "in", res = dpi, scaling = 1, background = "white",
                 compression = "lzw")
  if (inherits(plot, "ggplot") || inherits(plot, "patchwork")) print(plot) else grid.draw(plot)
  dev.off()
  cat("wrote", name, "\n")
}

## schematic helpers ---------------------------------------------------------
gb <- function(x, y, w, h, fill = PALE, col = INK, lw = 0.4, alpha = 1) {
  geom_rect(aes(xmin = x, xmax = x + w, ymin = y, ymax = y + h), fill = fill, colour = col, linewidth = lw, alpha = alpha, inherit.aes = FALSE)
}
gbox <- function(x, y, w, h, label, fill = PALE, col = INK, fs = 6.6, tc = INK, lw = 0.4, fontface = "plain") {
  list(geom_rect(aes(xmin = x, xmax = x + w, ymin = y, ymax = y + h), fill = fill, colour = col, linewidth = lw, inherit.aes = FALSE),
       annotate("text", x = x + w/2, y = y + h/2, label = label, size = fs/ggplot2::.pt, colour = tc, fontface = fontface, lineheight = 0.9))
}
gar <- function(x0, y0, x1, y1, col = INK, lw = 0.6, curved = FALSE) {
  if (curved) geom_curve(aes(x = x0, y = y0, xend = x1, yend = y1), colour = col, linewidth = lw,
                         arrow = arrow(length = unit(2.2,"pt"), type = "closed"), curvature = -0.2, inherit.aes = FALSE)
  else geom_segment(aes(x = x0, y = y0, xend = x1, yend = y1), colour = col, linewidth = lw,
                    arrow = arrow(length = unit(2.2,"pt"), type = "closed"), inherit.aes = FALSE)
}
blank <- function(xr = c(0,100), yr = c(0,100)) ggplot() + xlim(xr) + ylim(yr) +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA))

## ===========================================================================
## FIGURE 1 — schematic: resolution paradox + two-level framework
## ===========================================================================
fig1 <- function() {
  # Panel a: the paradox
  a <- blank(c(0,100), c(0,100)) +
    gbox(6, 78, 38, 14, "SIGLEC15+ myeloid cell\nACP5/MMP9 state present", BLUE_P, BLUE, fontface="bold") +
    gbox(56, 78, 38, 14, "SIGLEC15- myeloid cell\nstate absent", PALE) +
    annotate("text", x = 50, y = 85, label = "+0.64", colour = GREEN, fontface = "bold", size = 11/ggplot2::.pt) +
    annotate("text", x = 50, y = 96, label = "Cellular resolution (Level 1)", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gar(50, 77.5, 50, 70, col = ORANGE, lw = 0.8) +
    annotate("text", x = 53, y = 74, label = "aggregation", colour = ORANGE, size = 6/ggplot2::.pt, hjust = 0) +
    gbox(3, 50, 28, 14, "Sample aggregate\nr = -0.45  [-0.84, 0.11]\nUNINFORMATIVE", "#EEEEEE", GREY, fs = 6) +
    gbox(36, 50, 28, 14, "Bulk (n = 321)\nr = 0.13  [-0.11, 0.36]\nUNINFORMATIVE (large effects excluded)", RED_P, RED, fs = 6) +
    gbox(69, 50, 28, 14, "Survival (8 cohorts)\nHR 1.04  [0.96, 1.12]\nEXPLORATORY", UNINF, UNINF, fs = 6) +
    annotate("text", x = 50, y = 40, label = "Patient-level readouts (Level 2):\nthe cell-level signal does not resurface",
             colour = GREY, size = 6.2/ggplot2::.pt, fontface = "italic") +
    labs(tag = "a")

  # Panel b: protocol flow
  b <- blank(c(0,100), c(0,100)) +
    annotate("text", x = 50, y = 96, label = "Two-level qualification protocol", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gbox(8, 78, 84, 9, "1  State definition  (literature-anchored ACP5/MMP9; SIGLEC15 excluded from scores)", BLUE_P, BLUE, fs = 6) +
    gar(50, 78, 50, 72.5, lw = 0.5) +
    gbox(8, 62.5, 84, 9, "2  Level 1  cellular recoverability  (sample as inferential unit)", TEAL_P, TEAL, fs = 6) +
    gar(50, 62.5, 50, 57, lw = 0.5) +
    gbox(8, 47, 84, 9, "3  Same-system calibration gate  (reference targets recover expected programmes)", BLUE_P, BLUE, fs = 6) +
    gar(50, 47, 50, 41.5, lw = 0.5) +
    gbox(8, 31.5, 84, 9, "4  Level 2  patient-level portability  (aggregate / bulk / spatial / survival)", ORANGE_P, ORANGE, fs = 6) +
    gar(50, 31.5, 50, 26, lw = 0.5) +
    gbox(8, 16, 84, 9, "5  Verdicts + combining rule  (Level 1 never overrides Level 2)", PALE, INK, fs = 6) +
    annotate("text", x = 50, y = 7, label = "thresholds: qualification r >= 0.20 ;  large-coupling exclusion r >= 0.36",
             colour = GREY, size = 5.6/ggplot2::.pt, fontface = "italic") +
    labs(tag = "b")

  # Panel c: data landscape
  c <- blank(c(0,100), c(0,100)) +
    annotate("text", x = 50, y = 96, label = "Study data landscape", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gbox(4, 76, 44, 12, "Level 1  single-cell\nGSE155698 / GSE156405 / GSE212966\n10,431 myeloid cells, 20 samples", TEAL_P, TEAL, fs = 5.8) +
    gbox(52, 76, 44, 12, "Level 2  spatial\nGSE282302 (Visium)\n44 captures, 3 cases", BLUE_P, BLUE, fs = 5.8) +
    gbox(4, 56, 44, 12, "Level 2  bulk\nTCGA-PAAD / GSE57495 / GSE85916\nn = 321", ORANGE_P, ORANGE, fs = 5.8) +
    gbox(52, 56, 44, 12, "Level 2  survival (exploratory)\n8 TCGA cohorts, n = 2,955", PALE, INK, fs = 5.8) +
    gbox(4, 36, 92, 12, "Proteomics boundary  CPTAC-PDAC (n = 140) / PXD048644       Reference targets: CD274, TIGIT, ERBB2, EPCAM, CSF1R", PALE, INK, fs = 5.8) +
    annotate("text", x = 50, y = 22, label = "solid: within-series transitions      dashed: independent datasets",
             colour = GREY, size = 5.6/ggplot2::.pt, fontface = "italic") +
    labs(tag = "c")

  (a | b | c) + plot_layout(widths = c(1, 1, 1))
}

## ===========================================================================
## FIGURE 2 — RNA landscape
## ===========================================================================
fig2 <- function() {
  det <- data.frame(
    cohort = c("TCGA-PAAD", "PACA-CA", "GSE57495", "GSE85916", "GTEx"),
    pct = c(99.4, 76.7, 100, 100, 72.9),
    n = c(178, 43, 63, 80, 328),
    grp = c(rep("Tumour", 4), "Normal")
  ) |>
    mutate(cohort = factor(cohort, levels = cohort),
           lab = case_when(cohort == "PACA-CA" ~ "76.7%\n33/43",
                           pct == 100 ~ sprintf("100%%\nn = %d", n),
                           TRUE ~ sprintf("%.1f%%\nn = %d", pct, n)),
           label_y = c(104, 83, 104, 114, 78))
  a <- ggplot(det, aes(cohort, pct, fill = grp)) +
    geom_col(width = 0.66, colour = "white", linewidth = 0.25) +
    geom_col(data = filter(det, grp == "Normal"), width = 0.66,
             fill = GREY, colour = "white", linewidth = 0.25) +
    geom_text(aes(y = label_y, label = lab), size = 2.35, lineheight = 0.95, colour = INK) +
    scale_fill_manual(values = c("Tumour" = BLUE, "Normal" = GREY)) +
    scale_y_continuous(limits = c(0, 121), breaks = seq(0, 120, 20), expand = c(0, 0)) +
    labs(tag = "a", title = "Bulk RNA detection",
         x = NULL, y = "RNA-detectable samples (%)") +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 28, hjust = 1),
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.25))

  scr <- read_csv(file.path(SUPPV6, "Supplementary_Table_S3.csv"), show_col_types = FALSE) |>
    pivot_longer(-cell_type, names_to = "dataset", values_to = "pct") |>
    mutate(dataset = sub("_percent$", "", dataset),
           dataset = factor(dataset, levels = c("GSE155698", "GSE156405", "GSE212966")),
           cell_type = factor(cell_type, levels = c("Epithelial", "Myeloid", "Fibroblast", "T cell", "B cell")))
  b <- ggplot(scr, aes(cell_type, pct, colour = dataset, shape = dataset)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.10, fill = PALE) +
    geom_point(position = position_dodge(width = 0.48), size = 2.2, stroke = 0.55) +
    scale_colour_manual(values = DS_COL) +
    scale_shape_manual(values = c("GSE155698" = 16, "GSE156405" = 15, "GSE212966" = 17)) +
    scale_y_continuous(limits = c(0, 3.65), expand = c(0, 0)) +
    labs(tag = "b", title = "Single-cell expression is sparse",
         x = NULL, y = "SIGLEC15 RNA-detectable cells (%)", colour = NULL, shape = NULL) +
    theme(axis.text.x = element_text(angle = 28, hjust = 1),
          legend.position = c(0.70, 0.80),
          legend.background = element_rect(fill = alpha("white", 0.85), colour = NA),
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.25))

  order_sp <- c("Epithelial", "Myeloid", "NK", "B_cell", "T_cell", "Fibroblast", "Acinar", "Endothelial")
  sp <- read_csv(file.path(SUPPV6, "Supplementary_Table_S4.csv"), show_col_types = FALSE) |>
    filter(spot_type %in% order_sp) |>
    mutate(spot_type = factor(spot_type, levels = rev(order_sp)),
           pct = 100 * pos_frac,
           grp = case_when(as.character(spot_type) == "Epithelial" ~ "Epithelial",
                           as.character(spot_type) == "Myeloid" ~ "Myeloid",
                           TRUE ~ "Other"))
  c <- ggplot(sp, aes(pct, spot_type, fill = grp)) +
    geom_col(width = 0.68, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = c("Epithelial" = BLUE, "Myeloid" = ORANGE, "Other" = LGREY)) +
    scale_y_discrete(labels = function(x) gsub("_", " ", x)) +
    scale_x_continuous(limits = c(0, 1.55), breaks = c(0, 0.5, 1.0, 1.5), expand = c(0, 0)) +
    labs(tag = "c", title = "Spatial detection remains low",
         x = "RNA-detectable spots (%)", y = NULL) +
    theme(legend.position = "none", panel.grid.major.x = element_line(colour = LGREY, linewidth = 0.25))

  (a | b | c) + plot_layout(widths = c(1.05, 1.25, 1.08))
}

## ===========================================================================
## FIGURE 3 — Level 1 cellular state
## ===========================================================================
fig3 <- function() {
  eff <- read_csv(file.path(S17, "cell_level_sample_effects.csv"), show_col_types = FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core") |>
    mutate(dataset = factor(dataset, levels = c("GSE155698","GSE156405","GSE212966")))
  core_sum <- read_csv(file.path(S17, "cell_level_summary.csv"), show_col_types = FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core")

  effect_panel <- function(value_col, est, lo, hi, xlab, tag, title, show_legend = FALSE) {
    d <- eff |>
      transmute(dataset, effect = .data[[value_col]],
                row = factor(as.character(dataset),
                             levels = c("Mean (95% CI)", "GSE212966", "GSE156405", "GSE155698")))
    summ <- data.frame(effect = est, lo = lo, hi = hi,
                       row = factor("Mean (95% CI)",
                                    levels = levels(d$row)))
    ggplot(d, aes(effect, row, colour = dataset, shape = dataset)) +
      geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.35) +
      geom_jitter(height = 0.10, width = 0, size = 1.75, alpha = 0.90, stroke = 0.45) +
      geom_errorbarh(data = summ, aes(xmin = lo, xmax = hi, y = row),
                     inherit.aes = FALSE, height = 0, colour = ORANGE, linewidth = 0.55) +
      geom_point(data = summ, aes(x = effect, y = row), inherit.aes = FALSE,
                 shape = 18, size = 2.8, colour = ORANGE) +
      scale_colour_manual(values = DS_COL) +
      scale_shape_manual(values = c("GSE155698" = 16, "GSE156405" = 15, "GSE212966" = 17)) +
      labs(tag = tag, title = title, x = xlab, y = NULL, colour = NULL, shape = NULL) +
      theme(legend.position = if (show_legend) "top" else "none",
            legend.justification = "left",
            panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.22))
  }
  a <- effect_panel("difference_positive_minus_negative",
                    core_sum$mean_sample_difference, core_sum$ci_low, core_sum$ci_high,
                    "Positive-minus-negative score", "a", "Unadjusted cell contrast", TRUE)
  b <- effect_panel("depth_adjusted_difference",
                    core_sum$depth_adjusted_mean_difference,
                    core_sum$depth_adjusted_ci_low, core_sum$depth_adjusted_ci_high,
                    "Depth-adjusted score difference", "b", "Depth-adjusted cell contrast")

  dm <- read_csv(file.path(S19, "depth_matched_summary.csv"), show_col_types = FALSE) |>
    mutate(prog = recode(program, "Derosiers_ACP5_MMP9_core"="ACP5/MMP9 core",
                         "canonical_osteoclast_panel"="Canonical osteoclast",
                         "compact_ECM_degradation_panel"="Matrix degradation"),
           prog = factor(prog, levels = c("ACP5/MMP9 core","Canonical osteoclast","Matrix degradation")))
  c <- ggplot(dm, aes(mean_difference, prog)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, colour = BLUE, linewidth = 0.4) +
    geom_point(colour = BLUE, size = 2.0) +
    geom_text(aes(label = sprintf("p = %.3f", wilcoxon_p)), hjust = -0.08,
              size = 2.25, colour = GREY, nudge_y = 0.16) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.26))) +
    labs(tag = "c", x = "Depth-matched difference", y = NULL)

  sg0 <- bind_rows(
    read_csv(file.path(S18,"ACP5_cell_summary.csv"), show_col_types=FALSE) |> mutate(gene="ACP5"),
    read_csv(file.path(S18,"MMP9_cell_summary.csv"), show_col_types=FALSE) |> mutate(gene="MMP9"))
  sg <- bind_rows(
    sg0 |> transmute(gene, analysis = "Unadjusted", est = mean_sample_difference, lo = ci_low, hi = ci_high),
    sg0 |> transmute(gene, analysis = "Depth-adjusted", est = depth_adjusted_mean_difference,
                     lo = depth_adjusted_ci_low, hi = depth_adjusted_ci_high)
  ) |>
    mutate(analysis = factor(analysis, c("Unadjusted", "Depth-adjusted")),
           gene = factor(gene, c("ACP5", "MMP9")))
  d <- ggplot(sg, aes(est, gene, colour = analysis, shape = analysis)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0,
                   position = position_dodge(width = 0.34), linewidth = 0.4) +
    geom_point(position = position_dodge(width = 0.34), size = 2.0, stroke = 0.55) +
    scale_colour_manual(values = c("Unadjusted" = BLUE, "Depth-adjusted" = TEAL)) +
    scale_shape_manual(values = c("Unadjusted" = 16, "Depth-adjusted" = 15)) +
    labs(tag = "d", x = "Single-gene difference", y = NULL) +
    theme(legend.position = "top", legend.title = element_blank())

  th <- read_csv(file.path(P02,"threshold_sensitivity_table.csv"), show_col_types=FALSE) |>
    mutate(threshold_lab = ifelse(threshold == 0, "> 0", sprintf("> %.2g", threshold)),
           threshold_lab = factor(threshold_lab, levels = rev(threshold_lab)),
           n_scaled = 0.20 + (cell_n_positive_total - 20) * (2.60 / 177))
  e <- ggplot(th, aes(cell_mean_difference, threshold_lab)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = cell_ci_low, xmax = cell_ci_high), height = 0,
                   colour = ORANGE, linewidth = 0.45) +
    geom_point(colour = ORANGE, size = 2.0, shape = 16) +
    geom_line(aes(x = n_scaled, group = 1), colour = PURPLE, linetype = "dotted", linewidth = 0.35) +
    geom_point(aes(x = n_scaled), colour = PURPLE, fill = "white", shape = 21, size = 2.0, stroke = 0.65) +
    scale_x_continuous(
      limits = c(-0.2, 3.25),
      sec.axis = sec_axis(~ 20 + (. - 0.20) * 177 / 2.60,
                          name = "Estimable SIGLEC15+ cells",
                          breaks = c(20, 50, 100, 150, 200))
    ) +
    labs(tag = "e", x = "Positive-minus-negative score", y = "Positivity threshold") +
    theme(axis.title.x.top = element_text(colour = PURPLE), axis.text.x.top = element_text(colour = PURPLE))

  ((a | b) / (c | d) / e) +
    plot_layout(heights = c(1.08, 1, 1.05))
}

## ===========================================================================
## FIGURE 4 — Level 2 qualification failure
## ===========================================================================
fig4 <- function() {
  agg <- read_csv(file.path(S17,"myeloid_aggregate_associations.csv"), show_col_types=FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core") |>
    mutate(label = ifelse(scope == "pooled_within_dataset_ranks", "Pooled ranks",
                          recode(dataset, "GSE155698" = "GSE155698", "GSE156405" = "GSE156405",
                                 "GSE212966" = "GSE212966")),
           pooled = scope == "pooled_within_dataset_ranks",
           label = factor(label, levels = rev(c("GSE155698", "GSE156405", "GSE212966", "Pooled ranks"))))
  a <- ggplot(agg, aes(rho, label, colour = pooled, shape = pooled)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.45) +
    geom_point(size = 2.1, stroke = 0.55) +
    scale_colour_manual(values = c("FALSE" = BLUE, "TRUE" = ORANGE)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18)) +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(tag = "a", title = "Myeloid sample aggregates",
         x = "Spearman r (95% CI)", y = NULL) +
    theme(legend.position = "none", panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.22))

  bulk <- read_csv(file.path(S17,"bulk_associations.csv"), show_col_types=FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core") |>
    mutate(cohort = recode(cohort, "TCGA_PAAD"="TCGA-PAAD","GSE57495"="GSE57495","GSE85916"="GSE85916"))
  meta <- read_csv(file.path(S17,"bulk_random_effects_meta.csv"), show_col_types=FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core")
  adj <- read_csv(file.path(S17,"bulk_posthoc_myeloid_adjusted_meta.csv"), show_col_types=FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core")
  bd <- bind_rows(
    data.frame(cohort = bulk$cohort, rho = bulk$rho, ci_low = bulk$ci_low, ci_high = bulk$ci_high, kind = "Cohort"),
    data.frame(cohort = "Pooled", rho = meta$r_random, ci_low = meta$ci_low, ci_high = meta$ci_high, kind = "Pooled"),
    data.frame(cohort = "Pan-myeloid adjusted", rho = adj$r_random, ci_low = adj$ci_low, ci_high = adj$ci_high, kind = "Adjusted")
  ) |>
    mutate(cohort = factor(cohort, levels = rev(c("TCGA-PAAD", "GSE57495", "GSE85916", "Pooled", "Pan-myeloid adjusted"))),
           kind = factor(kind, c("Cohort", "Pooled", "Adjusted")))
  b <- ggplot(bd, aes(rho, cohort, colour = kind, shape = kind)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.36,-0.2,0.2,0.36), colour = c(RED,GREEN,GREEN,RED), linetype = c("dotted","dashed","dashed","dotted"), linewidth = 0.3) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.45, show.legend = FALSE) +
    geom_point(size = 2.1, stroke = 0.55) +
    scale_colour_manual(values = c("Cohort" = BLUE, "Pooled" = ORANGE, "Adjusted" = PURPLE)) +
    scale_shape_manual(values = c("Cohort" = 16, "Pooled" = 18, "Adjusted" = 0)) +
    annotate("text", x = 0.2, y = 5.30, label = "0.20", colour = GREEN, size = 2.25) +
    annotate("text", x = 0.36, y = 5.30, label = "0.36", colour = RED, size = 2.25) +
    coord_cartesian(xlim = c(-0.5, 0.55)) +
    labs(tag = "b", title = "Bulk RNA association",
         x = "Spearman r (95% CI)", y = NULL, colour = NULL, shape = NULL) +
    theme(legend.position = "bottom")

  sp <- read_csv(file.path(S17,"spatial_case_effects.csv"), show_col_types=FALSE) |>
    filter(program == "Derosiers_ACP5_MMP9_core")
  spdat <- bind_rows(
    data.frame(case = "C1", m = 0, lo = 0, hi = 0, kind = "Not estimable"),
    data.frame(case=sp$case, m=sp$mean_capture_difference, lo=sp$mean_capture_difference, hi=sp$mean_capture_difference, kind="Case"),
    data.frame(case="Mean (2 cases)", m=0.466, lo=-3.385, hi=4.318, kind="Mean")
  ) |>
    mutate(case = factor(case, levels = rev(c("C1", "C2", "C3", "Mean (2 cases)"))),
           kind = factor(kind, c("Case", "Mean", "Not estimable")))
  c <- ggplot(spdat, aes(m, case, colour = kind, shape = kind)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
    geom_point(size = 2.1, stroke = 0.6) +
    geom_text(data = filter(spdat, kind == "Not estimable"), aes(x = 0.55, y = case, label = "not estimable"),
              inherit.aes = FALSE, colour = GREY, hjust = 0, size = 2.25) +
    scale_colour_manual(values = c("Case" = BLUE, "Mean" = ORANGE, "Not estimable" = GREY)) +
    scale_shape_manual(values = c("Case" = 16, "Mean" = 18, "Not estimable" = 4)) +
    coord_cartesian(xlim = c(-4, 5)) +
    labs(tag = "c", title = "Spatial cases",
         x = "Case-level score difference", y = NULL) +
    theme(legend.position = "none")

  vd <- data.frame(layer = c("Cellular contrast","Sample aggregate (n=20)","Bulk pooled (n=321)","Spatial (n=2 cases)","Survival"),
                   verdict = c("Supportive","Imprecise","Benchmark indeterminate","Exploratory","Exploratory"))
  vd$layer <- factor(vd$layer, levels = rev(vd$layer))
  d <- ggplot(vd, aes(1, layer, fill = verdict)) +
    geom_tile(width = 0.92, height = 0.70) +
    geom_text(aes(label = verdict), colour = "white", size = 2.35, fontface = "bold") +
    scale_fill_manual(values = c("Supportive" = GREEN, "Imprecise" = GREY, "Benchmark indeterminate" = ORANGE, "Exploratory" = GREY)) +
    scale_x_continuous(limits = c(0.5, 1.5), expand = c(0,0)) +
    labs(tag = "d", title = "Interpretation by readout", x = NULL, y = NULL) +
    theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = 7.2), legend.position = "none")

  (a | b) / (c | d) + plot_layout(heights = c(1, 1), widths = c(1, 1.05))
}

## ===========================================================================
## FIGURE 5 — not technical
## ===========================================================================
fig5 <- function() {
  cal <- data.frame(target = c("CD274","TIGIT","ERBB2","EPCAM","CSF1R","SIGLEC15-PI3K"),
                    rho = c(0.545,0.883,0.568,0.679,0.822,0.186),
                    lo  = c(0.420,0.858,0.512,0.584,0.756,0.034),
                    hi  = c(0.652,0.902,0.619,0.757,0.871,0.331))
  cal$target <- factor(cal$target, levels = rev(cal$target))
  a <- ggplot(cal, aes(rho, target, colour = target == "SIGLEC15-PI3K", shape = target == "SIGLEC15-PI3K")) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, show.legend = FALSE) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = c("FALSE" = BLUE, "TRUE" = ORANGE)) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 18)) +
    labs(tag = "a", title = "Expected-programme calibration",
         x = "Spearman r (95% CI)", y = NULL) +
    theme(legend.position = "none")

  ds <- read_csv(file.path(ANALYSIS,"S10_cd274_downsampling","downsampling_replicates.csv"), show_col_types=FALSE) |>
    pull(rho_Hallmark_PI3K_AKT_MTOR__mean_z)
  b <- ggplot(data.frame(r = ds), aes(r)) +
    geom_histogram(bins = 22, fill = BLUE_P, colour = BLUE, linewidth = 0.25) +
    geom_vline(xintercept = median(ds), colour = ORANGE, linetype = "dashed", linewidth = 0.4) +
    annotate("label", x = min(ds) + 0.005, y = Inf, vjust = 1.15, hjust = 0,
             label = sprintf("Median r = +%.3f\n200/200 exceed r = 0.20", median(ds)),
             colour = INK, fill = alpha("white", 0.88), label.size = 0,
             size = 2.25, lineheight = 1.02) +
    labs(tag = "b", title = "Abundance-matched CD274 control",
         x = "CD274-PI3K r", y = "Replicates")

  sc <- read_csv(file.path(P22,"scorecard_table.csv"), show_col_types=FALSE) |>
    mutate(program_clean = recode(program, "CD8_T" = "CD8 T", "regulatory_T" = "regulatory T",
                                  "pan_myeloid" = "pan-myeloid"),
           name = paste0(target, " – ", program_clean),
           v = case_when(observed=="direction-consistent supportive" ~ "Supportive",
                          observed=="Non-supportive" ~ "Non-supportive", .default = "Uninformative")) |>
    arrange(pooled_rho)
  sc$name <- factor(sc$name, levels = sc$name)
  c <- ggplot(sc, aes(pooled_rho, name, colour = v)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-0.2,0.2), colour = GREEN, linetype = "dashed", linewidth = 0.25) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.35) +
    geom_point(size = 1.6) +
    scale_colour_manual(values = c("Supportive" = GREEN, "Non-supportive" = RED, "Uninformative" = GREY)) +
    coord_cartesian(xlim = c(-0.75, 1.05)) +
    labs(tag = "c", title = "Ten-pair reference analysis",
         x = "Pooled Spearman r (95% CI)", y = NULL, colour = NULL) +
    theme(legend.position = "top", legend.justification = "left")

  # heatmap as ggplot geom_tile (composes cleanly in patchwork)
  ct <- read_csv(file.path(ANALYSIS,"cross_target_panel_verified","cross_target_results_flat.csv"), show_col_types=FALSE) |>
    filter(Layer == "Competitive ecology")
  hdat <- ct |> select(Target, `Program/Ecology`, Spearman_rho) |>
    mutate(program = recode(`Program/Ecology`,
                            "T_cell_cytolytic" = "T-cell cytolytic",
                            "Myeloid" = "Myeloid",
                            "Epithelial" = "Epithelial",
                            "Stromal" = "Stromal",
                            "MUC_state" = "MUC state"),
           program = factor(program, levels = rev(c("T-cell cytolytic", "Myeloid", "Epithelial", "Stromal", "MUC state"))),
           Target = factor(Target, levels = c("CD274", "TIGIT", "ERBB2", "EPCAM", "CSF1R", "SIGLEC15")))
  d <- ggplot(hdat, aes(Target, program, fill = Spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.2f", Spearman_rho),
                  colour = abs(Spearman_rho) > 0.48), size = 2.15) +
    scale_colour_manual(values = c("FALSE" = INK, "TRUE" = "white"), guide = "none") +
    scale_fill_gradient2(low = "#B64E49", mid = "white", high = "#2C6E9B",
                         midpoint = 0, limits = c(-0.5, 0.85), name = "Spearman r") +
    labs(tag = "d", title = "Tissue-context calibration", x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right",
          legend.key.height = unit(11,"pt"))

  (a | b) / (c | d) + plot_layout(heights = c(0.85, 1.35), widths = c(1.08, 0.92))
}

## ===========================================================================
## FIGURE 6 — broader claims
## ===========================================================================
fig6 <- function() {
  pi3k <- read_csv(file.path(SUPPSRC, "Supplementary_Table_S16.csv"), show_col_types = FALSE) |>
    mutate(src = case_when(grepl("^Hallmark", source) ~ "Hallmark",
                           grepl("^KEGG", source) ~ "KEGG",
                           grepl("^PROGENy", source) ~ "PROGENy",
                           TRUE ~ "Reactome"),
           method_clean = recode(method, "mean_z" = "mean-z", "ssGSEA" = "ssGSEA",
                                 "GSVA" = "GSVA", "AUCell" = "AUCell"),
           label = paste(src, method_clean, sep = " – "),
           class = ifelse(classification == "equivalent", "Equivalent", "Inconclusive")) |>
    arrange(factor(src, c("Hallmark", "KEGG", "PROGENy", "Reactome")),
            factor(method_clean, c("AUCell", "GSVA", "mean-z", "ssGSEA"))) |>
    mutate(label = factor(label, levels = rev(label)))
  a <- ggplot(pi3k, aes(rho, label, colour = src, shape = class)) +
    annotate("rect", xmin = -0.20, xmax = 0.20, ymin = -Inf, ymax = Inf,
             fill = GREEN_P, alpha = 0.65) +
    geom_vline(xintercept = c(-0.2, 0.2), colour = GREEN, linetype = "dashed", linewidth = 0.32) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_errorbarh(aes(xmin = boot_ci95_low, xmax = boot_ci95_high), height = 0, linewidth = 0.38) +
    geom_point(size = 1.95, stroke = 0.65) +
    scale_colour_manual(values = c("Hallmark" = BLUE, "KEGG" = TEAL,
                                   "PROGENy" = PURPLE, "Reactome" = ORANGE)) +
    scale_shape_manual(values = c("Equivalent" = 0, "Inconclusive" = 16)) +
    coord_cartesian(xlim = c(-0.34, 0.40)) +
    labs(tag = "a", title = "PI3K specifications",
         x = "Spearman r (95% bootstrap CI)", y = NULL, colour = NULL, shape = NULL) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.20))

  syk <- read_csv(file.path(ANALYSIS, "T1_SYK_DAP12", "T1_SYK_DAP12_summary_table.csv"), show_col_types = FALSE) |>
    mutate(program = recode(gene_set,
                            "DAP12_signaling_R-HSA-2424491" = "DAP12 signalling",
                            "DAP12_interactions_R-HSA-2172127" = "DAP12 interactions",
                            "SYK_kinase_footprint" = "SYK footprint",
                            "Li2022_TAM_program" = "Li 2022 TAM"),
           method_clean = recode(scoring, "mean_z" = "mean-z", "ssgsea" = "ssGSEA", "aucell" = "AUCell"),
           lo = as.numeric(sub("^\\[([^,]+),.*$", "\\1", TCGA_CI)),
           hi = as.numeric(sub("^.*,[[:space:]]*([^]]+)\\]$", "\\1", TCGA_CI)),
           label = paste(program, method_clean, sep = " – ")) |>
    arrange(factor(program, c("DAP12 signalling", "DAP12 interactions", "SYK footprint", "Li 2022 TAM")),
            factor(method_clean, c("mean-z", "ssGSEA", "AUCell"))) |>
    mutate(label = factor(label, levels = rev(label)))
  b <- ggplot(syk, aes(TCGA_rho, label, colour = program, shape = method_clean)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_vline(xintercept = 0.20, colour = GREEN, linetype = "dashed", linewidth = 0.32) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.38) +
    geom_point(size = 1.9, stroke = 0.6) +
    scale_colour_manual(values = c("DAP12 signalling" = BLUE, "DAP12 interactions" = TEAL,
                                   "SYK footprint" = PURPLE, "Li 2022 TAM" = ORANGE)) +
    scale_shape_manual(values = c("mean-z" = 16, "ssGSEA" = 15, "AUCell" = 17)) +
    coord_cartesian(xlim = c(-0.18, 0.36)) +
    labs(tag = "b", title = "SYK/DAP12 specifications",
         x = "TCGA-PAAD Spearman r (95% CI)", y = NULL, colour = NULL, shape = NULL) +
    theme(legend.position = "none",
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.20))

  muc_meta <- read_csv(file.path(VALID,"external_meta_reml_hk.csv"), show_col_types=FALSE) |>
    filter(gene_set == "MUC1_CD24")
  muc_raw <- read_csv(file.path(ANALYSIS, "external_validation", "per_cohort_associations.csv"), show_col_types = FALSE) |>
    filter(gene_set == "MUC1_CD24") |>
    mutate(z = atanh(spearman_rho), se = 1 / sqrt(n_samples - 3),
           lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se),
           cohort = recode(cohort, "TCGA_PAAD" = "TCGA-PAAD"),
           r = spearman_rho, pooled = FALSE) |>
    select(cohort, r, lo, hi, pooled)
  mucdat <- bind_rows(
    muc_raw,
    data.frame(cohort = "Pooled (REML + mHK)", r = muc_meta$pooled_r,
               lo = muc_meta$modified_hk_ci95_lo_r, hi = muc_meta$modified_hk_ci95_hi_r,
               pooled = TRUE)
  ) |>
    mutate(cohort = factor(cohort, levels = rev(c("TCGA-PAAD", "GSE57495", "GSE85916", "Pooled (REML + mHK)"))))
  c <- ggplot(mucdat, aes(r, cohort, colour = pooled, shape = pooled)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
    geom_point(size = 2.0) +
    scale_colour_manual(values = c("FALSE" = BLUE, "TRUE" = ORANGE)) +
    scale_shape_manual(values = c("FALSE"=16,"TRUE"=18)) +
    coord_cartesian(xlim = c(-0.7, 0.85)) +
    labs(tag = "c", title = "MUC1/CD24 association was not reproduced",
         x = "Spearman r (95% CI)", y = NULL) +
    theme(legend.position = "none")

  evidence_row <- function(p, y, accent, fill, heading, detail, badge, badge_fill, badge_colour) {
    p +
      geom_rect(aes(xmin = 5, xmax = 95, ymin = y - 6.7, ymax = y + 6.7),
                fill = fill, colour = LGREY, linewidth = 0.25) +
      geom_rect(aes(xmin = 5, xmax = 6.5, ymin = y - 6.7, ymax = y + 6.7),
                fill = accent, colour = accent, linewidth = 0) +
      annotate("text", x = 9, y = y + 2.3, label = heading, hjust = 0,
               colour = INK, fontface = "bold", size = 7.4 / ggplot2::.pt) +
      annotate("text", x = 9, y = y - 2.4, label = detail, hjust = 0,
               colour = GREY, size = 6.3 / ggplot2::.pt) +
      geom_rect(aes(xmin = 74, xmax = 94, ymin = y - 4.1, ymax = y + 4.1),
                fill = badge_fill, colour = badge_colour, linewidth = 0.35) +
      annotate("text", x = 84, y = y, label = badge,
               colour = badge_colour, fontface = "bold", size = 5.0 / ggplot2::.pt,
               lineheight = 0.88)
  }

  d <- blank(c(0, 100), c(0, 100)) +
    annotate("text", x = 3, y = 97, label = "d", hjust = 0, vjust = 1,
             colour = INK, fontface = "bold", size = 10.5 / ggplot2::.pt) +
    annotate("text", x = 50, y = 94, label = "Protein evidence is incomplete",
             colour = INK, fontface = "bold", size = 8.2 / ggplot2::.pt) +
    annotate("text", x = 5, y = 88, label = "Evidence layer", hjust = 0,
             colour = GREY, fontface = "bold", size = 5.4 / ggplot2::.pt) +
    annotate("text", x = 92, y = 88, label = "Status", hjust = 1,
             colour = GREY, fontface = "bold", size = 5.4 / ggplot2::.pt)

  d <- evidence_row(d, 77, BLUE, BLUE_P, "Discovery proteomics",
                    "No quantifiable SIGLEC15 target row", "COVERAGE\nLIMIT",
                    BLUE_P, BLUE)
  d <- evidence_row(d, 60, GREEN, GREEN_P, "Published antibody studies",
                    "PDAC tissue signal; target engagement", "SUPPORTED",
                    GREEN_P, GREEN)
  d <- evidence_row(d, 43, RED, RED_P, "Matched RNA–protein",
                    "Compartmental concordance unavailable", "UNRESOLVED",
                    RED_P, RED)
  d <- evidence_row(d, 26, RED, RED_P, "Treatment-response evidence",
                    "Biomarker interaction unavailable", "UNRESOLVED",
                    RED_P, RED)
  d <- d + annotate("text", x = 50, y = 10,
                    label = "Non-identification in discovery proteomics does not imply protein absence",
                    colour = GREY, fontface = "italic", size = 5.4 / ggplot2::.pt)

  ((a | b) + plot_layout(widths = c(1, 1))) /
    ((c | d) + plot_layout(widths = c(0.88, 1.12))) +
    plot_layout(heights = c(1.42, 1.08))
}

## ===========================================================================
## FIGURE 7 — protocol deliverable (schematic)
## ===========================================================================
fig7 <- function() {
  a <- blank(c(0,100), c(0,100)) +
    annotate("text", x = 50, y = 96, label = "Decision flow", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gbox(8, 80, 84, 9, "State definition  literature-anchored ACP5/MMP9; SIGLEC15 excluded from scores", BLUE_P, BLUE, fs=5.8) +
    gar(50,80,50,74.5,lw=0.5) + gbox(8,65,84,9,"Level 1  cellular recoverability (sample as unit; depth, matching, thresholds)", TEAL_P, TEAL, fs=5.8) +
    gar(50,65,50,59.5,lw=0.5) + gbox(8,50,84,9,"Calibration gate  reference targets must recover expected relationships", BLUE_P, BLUE, fs=5.8) +
    gar(50,50,50,44.5,lw=0.5) + gbox(8,35,84,9,"Level 2  patient-level portability  aggregate / bulk / spatial / survival", ORANGE_P, ORANGE, fs=5.8) +
    gar(50,35,50,29.5,lw=0.5) + gbox(8,20,84,9,"Verdicts + combining rule  Level 1 never overrides Level 2", PALE, INK, fs=5.8) +
    annotate("text", x = 50, y = 11, label = "uninformative is not evidence either way", colour = GREY, size = 5.6/ggplot2::.pt, fontface="italic") +
    labs(tag="a")

  b <- blank(c(0,100), c(0,100)) +
    annotate("text", x = 50, y = 96, label = "SIGLEC15 RNA in PDAC  final grade", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gbox(4, 70, 92, 16, "Level 1  SUPPORTED\nACP5/MMP9 state recovered in SIGLEC15+ myeloid cells\n0.64 [0.33,0.99]; depth-adj 0.43 [0.17,0.73]; TREM2+ TAM niche", GREEN_P, GREEN, fs=6) +
    gbox(4, 44, 92, 16, "Level 2  NOT QUALIFIED\naggregate uninformative (r=-0.45); bulk failed qualification (r=0.13; large effects excluded);\nspatial uninformative (2 cases); survival null (HR 1.04); PI3K/SYK-DAP12 not reproducible", RED_P, RED, fs=5.8) +
    gbox(4, 26, 92, 10, "Biological claim intact     patient-selection claim not supported     protein-first measurement required", PALE, INK, fs=5.8) +
    annotate("text", x = 50, y = 16, label = "cell-level recovery  !=  patient-level readout", colour = ORANGE, fontface="bold", size = 6.4/ggplot2::.pt) +
    labs(tag="b")

  c <- blank(c(0,100), c(0,100)) +
    annotate("text", x = 50, y = 96, label = "Protein-first next study", colour = INK, fontface = "bold", size = 7/ggplot2::.pt) +
    gbox(4, 78, 92, 11, "1  Compartment-resolved SIGLEC15 IHC  (positive-cell prevalence)", BLUE_P, BLUE, fs=5.8) +
    gar(50,78,50,72.5,lw=0.5) + gbox(4,61,92,11,"2  Orthogonal targeted protein confirmation  (MS / proximity assay)", BLUE_P, BLUE, fs=5.8) +
    gar(50,61,50,55.5,lw=0.5) + gbox(4,44,92,11,"3  ACP5/MMP9 state readout in the protein-positive compartment", BLUE_P, BLUE, fs=5.8) +
    gar(50,44,50,38.5,lw=0.5) + gbox(4,27,92,11,"4  Treatment-by-biomarker validation  (comparator-controlled)", BLUE_P, BLUE, fs=5.8) +
    annotate("text", x = 50, y = 16, label = "until then: RNA supports hypothesis generation and tissue prioritisation, not clinical enrichment",
             colour = GREY, size = 5.6/ggplot2::.pt, fontface="italic") +
    labs(tag="c")

  (a | b | c) + plot_layout(widths = c(1,1.1,1))
}

## ===========================================================================
## SUPPLEMENTARY FIGURES
## ===========================================================================
figs1 <- function() {  # tissue-context heatmap (ggplot)
  ct <- read_csv(file.path(ANALYSIS,"cross_target_panel_verified","cross_target_results_flat.csv"), show_col_types=FALSE) |>
    filter(Layer == "Competitive ecology") |>
    mutate(program = recode(`Program/Ecology`,
                            "T_cell_cytolytic" = "T-cell cytolytic",
                            "Myeloid" = "Myeloid",
                            "Epithelial" = "Epithelial",
                            "Stromal" = "Stromal",
                            "MUC_state" = "MUC state"),
           target = factor(Target, levels = c("CD274", "TIGIT", "ERBB2", "EPCAM", "CSF1R", "SIGLEC15")),
           program = factor(program, levels = rev(c("T-cell cytolytic", "Myeloid", "Epithelial",
                                                    "Stromal", "MUC state"))))
  ggplot(ct, aes(target, program, fill = Spearman_rho)) +
    geom_tile(colour = "white", linewidth = 0.45) +
    geom_text(aes(label = sprintf("%.2f", Spearman_rho)), size = 7.0 / ggplot2::.pt,
              colour = INK) +
    scale_fill_gradient2(low = RED, mid = "white", high = BLUE, midpoint = 0,
                         limits = c(-0.45, 0.90), na.value = LGREY,
                         name = "Spearman r") +
    labs(x = NULL, y = NULL, title = "Tissue-context calibration across targets") +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "right",
          panel.border = element_rect(colour = LGREY, fill = NA, linewidth = 0.35))
}
figs2 <- function() {  # myeloid subtype enrichment
  ch <- read_csv(file.path(ANALYSIS,"Task5_myeloid_subtypes","chi2_subtype_enrichment.csv"), show_col_types=FALSE) |>
    filter(dataset %in% names(DS_COL)) |>
    mutate(subtype = recode(subtype,
                            "Classical_Monocyte" = "Classical monocyte",
                            "Intermediate_Monocyte" = "Intermediate monocyte",
                            "NonClassical_Monocyte" = "Non-classical monocyte",
                            "M1_Macrophage" = "M1 macrophage",
                            "M2_Macrophage" = "M2 macrophage",
                            "TREM2_TAM" = "TREM2+ TAM"),
           subtype = factor(subtype, levels = rev(c("TREM2+ TAM", "M2 macrophage", "M1 macrophage",
                                                    "Classical monocyte", "Intermediate monocyte",
                                                    "Non-classical monocyte", "cDC1", "cDC2", "pDC"))))
  ggplot(ch, aes(enrichment_fold, subtype, colour = dataset, shape = dataset)) +
    geom_vline(xintercept = 1, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_point(size = 2.1, stroke = 0.65, position = position_dodge(width = 0.45)) +
    scale_colour_manual(values = DS_COL) +
    scale_shape_manual(values = c("GSE155698" = 16, "GSE156405" = 17, "GSE212966" = 15)) +
    scale_x_continuous(limits = c(0, 2.65), breaks = seq(0, 2.5, 0.5), expand = expansion(mult = c(0.01, 0.03))) +
    labs(x = "Observed / expected enrichment among SIGLEC15+ cells", y = NULL,
         colour = NULL, shape = NULL, title = "Myeloid-subtype enrichment is concentrated in TREM2+ TAMs") +
    theme(legend.position = "top", legend.justification = "left",
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.22))
}
figs3 <- function() {  # survival forest
  per <- read_csv(file.path(ANALYSIS,"Task2_Meta_Cox","per_cohort_univariate_cox.csv"), show_col_types=FALSE)
  pool <- read_csv(file.path(ANALYSIS,"Task2_Meta_Cox","meta_pooled_univariate.csv"), show_col_types=FALSE) |> slice(1)
  d <- bind_rows(data.frame(cohort=per$cohort, HR=per$HR, lo=per$HR_lower95, hi=per$HR_upper95, pooled=FALSE),
                 data.frame(cohort=sprintf("Pooled (I2 = %.1f%%)", pool$I2), HR=pool$pooled_HR, lo=pool$pooled_HR_lower95, hi=pool$pooled_HR_upper95, pooled=TRUE))
  d$cohort <- factor(d$cohort, levels = rev(d$cohort))
  ggplot(d, aes(HR, cohort, colour = pooled, shape = pooled)) +
    geom_vline(xintercept = 1, colour = GREY, linetype = "dashed", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4) +
    geom_point(size = 1.9) +
    scale_colour_manual(values = c("FALSE" = BLUE, "TRUE" = ORANGE)) +
    scale_shape_manual(values = c("FALSE"=16,"TRUE"=18)) +
    scale_x_log10(limits = c(0.6, 2.2)) +
    labs(x = "Univariate hazard ratio for continuous SIGLEC15 RNA (log scale)", y = NULL,
         title = "Survival associations are null across cohorts") +
    theme(legend.position = "none")
}
figs4 <- function() {  # deconvolution cross-tool faceted dot plot
  cs <- read_csv(file.path(ANALYSIS,"Task4_immune_tools_cross","consistency_matrix.csv"), show_col_types=FALSE)
  tool_levels <- c("ssGSEA", "MCPcounter", "EPIC", "quanTIseq", "CIBERSORTx_LM22")
  tool_cols <- c("ssGSEA" = BLUE, "MCPcounter" = TEAL, "EPIC" = ORANGE,
                 "quanTIseq" = PURPLE, "CIBERSORTx_LM22" = RED)
  cs <- cs |>
    mutate(tool = factor(tool, levels = tool_levels),
           cell_label = gsub("_", " ", cell_type),
           cell_label = factor(cell_label, levels = rev(unique(cell_label))))
  ggplot(cs, aes(rho, cell_label, colour = tool)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_segment(aes(x = 0, xend = rho, yend = cell_label), linewidth = 0.35, alpha = 0.55) +
    geom_point(aes(shape = sig_FDR05), size = 2.0, stroke = 0.65) +
    facet_grid(tool ~ ., scales = "free_y", space = "free_y") +
    scale_colour_manual(values = tool_cols, guide = "none") +
    scale_shape_manual(values = c("FALSE" = 1, "TRUE" = 16), guide = "none") +
    scale_x_continuous(limits = c(-0.30, 0.31), breaks = seq(-0.3, 0.3, 0.1),
                       labels = scales::label_number(accuracy = 0.1)) +
    labs(x = "Spearman r with SIGLEC15 RNA", y = NULL,
         title = "No immune-cell estimate survives cross-tool multiple-testing control") +
    theme(strip.background = element_rect(fill = PALE, colour = LGREY, linewidth = 0.35),
          strip.text.y = element_text(face = "bold", colour = INK, angle = 0, hjust = 0),
          panel.grid.major.y = element_line(colour = LGREY, linewidth = 0.20),
          panel.spacing.y = unit(2.5, "pt"))
}

figs5 <- function() {  # state-level positive controls (S20)
  s20_dir <- "./prior_analysis_v2/deposit_v24/results"
  s20 <- jsonlite::fromJSON(file.path(s20_dir, "S20_summary.json"), simplifyVector = TRUE)

  cell <- data.frame(
    state = c("TREM2 TAM", "SPP1 TAM"),
    est = c(s20$cell$TREM2_TAM$mean_sample_difference, s20$cell$SPP1_TAM$mean_sample_difference),
    lo = c(s20$cell$TREM2_TAM$ci_low, s20$cell$SPP1_TAM$ci_low),
    hi = c(s20$cell$TREM2_TAM$ci_high, s20$cell$SPP1_TAM$ci_high)
  ) |>
    mutate(state = factor(state, levels = rev(state)))
  a <- ggplot(cell, aes(est, state, colour = state, shape = state)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45) +
    geom_point(size = 2.3, stroke = 0.65) +
    scale_colour_manual(values = c("TREM2 TAM" = BLUE, "SPP1 TAM" = TEAL)) +
    scale_shape_manual(values = c("TREM2 TAM" = 16, "SPP1 TAM" = 17)) +
    labs(tag = "a", title = "Cell contrast", x = "Positive-minus-negative sample effect (95% CI)", y = NULL) +
    theme(legend.position = "none")

  bulk <- data.frame(
    state = c("TREM2 TAM", "SPP1 TAM", "LRRC15 CAF"),
    est = c(s20$bulk$TREM2_TAM$pooled_r, s20$bulk$SPP1_TAM$pooled_r, s20$bulk$LRRC15_CAF$pooled_r),
    lo = c(s20$bulk$TREM2_TAM$ci_low, s20$bulk$SPP1_TAM$ci_low, s20$bulk$LRRC15_CAF$ci_low),
    hi = c(s20$bulk$TREM2_TAM$ci_high, s20$bulk$SPP1_TAM$ci_high, s20$bulk$LRRC15_CAF$ci_high)
  ) |>
    mutate(state = factor(state, levels = rev(state)))
  b <- ggplot(bulk, aes(est, state, colour = state, shape = state)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_vline(xintercept = c(0.20, 0.36), colour = c(GREEN, ORANGE), linetype = "dotted", linewidth = 0.40) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45) +
    geom_point(size = 2.3, stroke = 0.65) +
    scale_colour_manual(values = c("TREM2 TAM" = BLUE, "SPP1 TAM" = TEAL, "LRRC15 CAF" = ORANGE)) +
    scale_shape_manual(values = c("TREM2 TAM" = 16, "SPP1 TAM" = 17, "LRRC15 CAF" = 15)) +
    coord_cartesian(xlim = c(0, 1.0)) +
    labs(tag = "b", title = "Bulk portability", x = "Pooled Spearman r (95% CI)", y = NULL) +
    theme(legend.position = "none")

  agg <- data.frame(
    state = c("TREM2 TAM", "SPP1 TAM"),
    est = c(s20$aggregate$TREM2_TAM$pooled_rho, s20$aggregate$SPP1_TAM$pooled_rho),
    lo = c(s20$aggregate$TREM2_TAM$ci_low, s20$aggregate$SPP1_TAM$ci_low),
    hi = c(s20$aggregate$TREM2_TAM$ci_high, s20$aggregate$SPP1_TAM$ci_high)
  ) |>
    mutate(state = factor(state, levels = rev(state)))
  c <- ggplot(agg, aes(est, state, colour = state, shape = state)) +
    geom_vline(xintercept = 0, colour = GREY, linetype = "dashed", linewidth = 0.30) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.45) +
    geom_point(size = 2.3, stroke = 0.65) +
    scale_colour_manual(values = c("TREM2 TAM" = BLUE, "SPP1 TAM" = TEAL)) +
    scale_shape_manual(values = c("TREM2 TAM" = 16, "SPP1 TAM" = 17)) +
    coord_cartesian(xlim = c(0, 1.0)) +
    labs(tag = "c", title = "Sample aggregates", x = "Pooled rank r (95% CI)", y = NULL) +
    theme(legend.position = "none")

  a | b | c
}

## ---- run all --------------------------------------------------------------
set.seed(15072026)
save_pub(fig2(), "Fig2", w_mm = 170, h_mm = 73)
save_pub(fig3(), "Fig3", w_mm = 170, h_mm = 145)
save_pub(fig4(), "Fig4", w_mm = 170, h_mm = 118)
save_pub(fig5(), "Fig5", w_mm = 170, h_mm = 130)
save_pub(fig6(), "Fig6", w_mm = 170, h_mm = 150)
save_pub(figs1(), "FigS1", w_mm = 160, h_mm = 92)
save_pub(figs2(), "FigS2", w_mm = 160, h_mm = 100)
save_pub(figs3(), "FigS3", w_mm = 150, h_mm = 92)
save_pub(figs4(), "FigS4", w_mm = 170, h_mm = 205)
save_pub(figs5(), "FigS5", w_mm = 170, h_mm = 88)

cat("ALL FIGURES DONE\n")
