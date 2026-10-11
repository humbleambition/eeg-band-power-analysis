# 05_figures.R ---------------------------------------------------------------
# All figures for the README and report, written to figures/*.png

source("R/utils.R")
suppressPackageStartupMessages({ library(patchwork); library(pROC) })
dir.create("figures", showWarnings = FALSE)

psd  <- readRDS(file.path("data", "derived", "psd.rds"))
bp   <- readr::read_csv(file.path("results", "band_power_long.csv"), show_col_types = FALSE) |>
  mutate(band = factor(band, levels(BANDS$band)), region = factor(region, REGION_LEVELS),
         group = factor(group, c("Control", "Schizophrenia")))
uni  <- readr::read_csv(file.path("results", "channel_band_stats.csv"), show_col_types = FALSE) |>
  mutate(band = factor(band, levels(BANDS$band)), channel = factor(channel, CHANNELS))
apf  <- readr::read_csv(file.path("results", "alpha_peak.csv"), show_col_types = FALSE)
apfs <- readr::read_csv(file.path("results", "alpha_peak_stats.csv"), show_col_types = FALSE)
n_ctl <- n_distinct(bp$subject[bp$group == "Control"])
n_scz <- n_distinct(bp$subject[bp$group == "Schizophrenia"])
grp_lab <- c(Control = sprintf("Control (n = %d)", n_ctl),
             Schizophrenia = sprintf("Schizophrenia spectrum (n = %d)", n_scz))
band_lines <- c(4, 8, 13)
band_labs  <- tibble(x = c(2.3, 6, 10.5, 21.5), lab = c("delta", "theta", "alpha", "beta"))

# ---- Fig 1: power spectra --------------------------------------------------
spec_summary <- function(chans) {
  psd |>
    filter(channel %in% chans, freq >= 1, freq <= 30) |>
    group_by(subject, group, freq) |>
    summarise(lp = mean(log10(psd)), .groups = "drop") |>
    group_by(group, freq) |>
    summarise(m = mean(lp), se = sd(lp) / sqrt(n()), .groups = "drop")
}
spec_plot <- function(d, title) {
  ggplot(d, aes(freq, m, colour = group, fill = group)) +
    geom_vline(xintercept = band_lines, linetype = "dashed", colour = scales::alpha(INK, 0.14)) +
    geom_ribbon(aes(ymin = m - 1.96 * se, ymax = m + 1.96 * se), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 0.8, lineend = "round") +
    geom_text(data = band_labs, aes(x = x, y = Inf, label = lab), inherit.aes = FALSE,
              vjust = 1.4, size = 3.3, colour = scales::alpha(INK, 0.5)) +
    scale_colour_manual(values = GROUP_COLS, labels = grp_lab) +
    scale_fill_manual(values = GROUP_COLS, labels = grp_lab) +
    scale_x_continuous(breaks = c(1, 4, 8, 13, 20, 30)) +
    labs(x = "Frequency (Hz)", y = expression(log[10]~"power density"), subtitle = title) +
    theme_portfolio()
}
p1 <- spec_plot(spec_summary(CHANNELS), "Mean of all 16 channels") +
      spec_plot(spec_summary(c("O1", "O2")), "Occipital (O1, O2)") +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Adolescents with schizophrenia-spectrum disorders show a smaller, slower alpha peak",
    subtitle = "Resting EEG, eyes closed. Group mean with 95% CI; Welch PSD, 4-s Hann windows",
    theme = theme_portfolio()) & theme(legend.position = "top")
save_fig(p1, "fig1_power_spectra.png", width = 11, height = 5)

# ---- Fig 2: effect-size heatmap --------------------------------------------
meas_lab <- c(abs_log10 = "Absolute power (log10)", rel = "Relative power (% of 1-30 Hz)")
hm <- uni |>
  mutate(measure = factor(meas_lab[measure], meas_lab),
         star = case_when(p_fdr < .001 ~ "***", p_fdr < .01 ~ "**", p_fdr < .05 ~ "*", TRUE ~ ""),
         channel = factor(channel, rev(CHANNELS)))
lim <- max(abs(hm$g))
p2 <- ggplot(hm, aes(band, channel, fill = g)) +
  geom_tile(colour = BG, linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%s", g, star)), size = 3, colour = INK, family = "mono") +
  facet_wrap(~measure) +
  scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]], high = DIVERGING[["high"]],
                       midpoint = 0, limits = c(-lim, lim), name = "Hedges' g") +
  labs(x = NULL, y = NULL,
       title = "Alpha is lower and slow-wave share is higher in schizophrenia",
       subtitle = "Hedges' g, schizophrenia minus control (positive = higher in schizophrenia). * BH-FDR q < .05, ** < .01, *** < .001 (64 tests per panel)") +
  theme_portfolio() +
  theme(panel.grid.major.y = element_blank(), legend.position = "right",
        legend.title = element_text(colour = scales::alpha(INK, 0.7), size = 10))
save_fig(p2, "fig2_effect_size_heatmap.png", width = 11, height = 7)

# ---- Fig 3: schematic scalp maps of effect size ----------------------------
head_circle <- tibble(t = seq(0, 2 * pi, length.out = 200), x = 1.12 * cos(t), y = 1.12 * sin(t))
nose <- tibble(x = c(-0.12, 0, 0.12), y = c(1.11, 1.25, 1.11))
topo <- uni |>
  mutate(measure = factor(meas_lab[measure], meas_lab)) |>
  left_join(ELEC_XY, by = "channel")
p3 <- ggplot(topo, aes(x, y)) +
  geom_path(data = head_circle, colour = scales::alpha(INK, 0.3)) +
  geom_path(data = nose, colour = scales::alpha(INK, 0.3)) +
  geom_point(aes(fill = g, size = abs(g), stroke = ifelse(p_fdr < .05, 1, 0.2)),
             shape = 21, colour = INK) +
  geom_text(aes(y = y + 0.19, label = channel), size = 2.3, colour = scales::alpha(INK, 0.55)) +
  facet_grid(measure ~ band, switch = "y") +
  scale_fill_gradient2(low = DIVERGING[["low"]], mid = DIVERGING[["mid"]], high = DIVERGING[["high"]],
                       midpoint = 0, limits = c(-lim, lim), name = "Hedges' g") +
  scale_size_continuous(range = c(1.5, 7), guide = "none") +
  coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.2, 1.35)) +
  labs(title = "Where the group differences sit on the scalp",
       subtitle = "Schematic 10-20 layout (nose up). Colour/size = Hedges' g; thick outline = FDR q < .05") +
  theme_void(base_size = 12) +
  theme(plot.background = element_rect(fill = BG, colour = NA),
        strip.text = element_text(face = "bold", colour = INK, margin = margin(4, 4, 4, 4)),
        strip.text.y.left = element_text(angle = 90),
        plot.title = element_text(face = "bold", colour = INK, size = 14),
        plot.subtitle = element_text(colour = scales::alpha(INK, 0.65), size = 11),
        legend.title = element_text(colour = scales::alpha(INK, 0.7), size = 10),
        plot.margin = margin(12, 12, 10, 12))
save_fig(p3, "fig3_scalp_effect_maps.png", width = 11, height = 6.2)

# ---- Fig 4: relative power by region ---------------------------------------
reg <- bp |>
  group_by(subject, group, region, band) |>
  summarise(rel = mean(rel) * 100, .groups = "drop")
p4 <- ggplot(reg, aes(region, rel, fill = group, colour = group)) +
  geom_boxplot(alpha = 0.35, outlier.shape = NA, width = 0.65,
               position = position_dodge(0.75), linewidth = 0.4) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75, seed = 1),
             size = 0.7, alpha = 0.55) +
  facet_wrap(~band, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = GROUP_COLS, labels = grp_lab) +
  scale_colour_manual(values = GROUP_COLS, labels = grp_lab) +
  scale_y_continuous(labels = function(v) paste0(v, "%"), limits = c(0, NA)) +
  labs(x = NULL, y = "Relative power (share of 1-30 Hz)",
       title = "Relative band power by scalp region",
       subtitle = "Each point is one adolescent (mean of the region's electrodes)") +
  theme_portfolio() +
  theme(axis.text.x = element_text(size = 10))
save_fig(p4, "fig4_relative_power_by_region.png", width = 11, height = 8)

# ---- Fig 5: alpha peak frequency -------------------------------------------
apf_all <- apf
apf <- mutate(apf, group = factor(group, c("Control", "Schizophrenia")))
p5 <- ggplot(filter(apf, !is.na(alpha_peak_hz)), aes(group, alpha_peak_hz, fill = group, colour = group)) +
  geom_violin(alpha = 0.25, colour = NA, width = 0.8) +
  geom_boxplot(width = 0.18, alpha = 0.6, outlier.shape = NA, colour = INK, linewidth = 0.35) +
  geom_point(position = position_jitter(width = 0.08, seed = 2), size = 1.3, alpha = 0.7) +
  scale_fill_manual(values = GROUP_COLS) + scale_colour_manual(values = GROUP_COLS) +
  scale_x_discrete(labels = c(Control = sprintf("Control (n = %d)", apfs$n_control),
                              Schizophrenia = sprintf("Schizophrenia spectrum (n = %d)", apfs$n_scz))) +
  labs(x = NULL, y = "Alpha peak frequency (Hz)",
       title = sprintf("Alpha peak is %.1f Hz slower in schizophrenia", apfs$mean_control - apfs$mean_scz),
       subtitle = sprintf("Posterior spectrum (O1, O2, P3, P4, Pz), 7-13 Hz; %d without a clear peak omitted\nWelch t(%.1f) = %.2f, p = %.4f; Hedges' g = %.2f [95%% CI %.2f, %.2f]",
                          sum(is.na(apf$alpha_peak_hz)), apfs$df, apfs$t, apfs$p, apfs$g, apfs$g_lo, apfs$g_hi)) +
  theme_portfolio() + theme(legend.position = "none")
save_fig(p5, "fig5_alpha_peak_frequency.png", width = 8.5, height = 5.2)

# ---- Fig 6: classifier -----------------------------------------------------
if (file.exists(file.path("results", "classifier_oof.csv"))) {
  oof  <- readr::read_csv(file.path("results", "classifier_oof.csv"), show_col_types = FALSE)
  perm <- readr::read_csv(file.path("results", "classifier_perm.csv"), show_col_types = FALSE)
  summ <- readr::read_csv(file.path("results", "classifier_summary.csv"), show_col_types = FALSE)
  roc_df <- oof |>
    pivot_longer(c(p_enet, p_baseline), names_to = "model", values_to = "p") |>
    group_by(repeat_id, model) |>
    group_modify(~{
      r <- pROC::roc(.x$y, .x$p, quiet = TRUE, direction = "<")
      tibble(fpr = 1 - r$specificities, tpr = r$sensitivities)
    }) |> ungroup() |>
    mutate(model = recode(model, p_enet = "Elastic net, 128 features",
                          p_baseline = "Baseline, global relative alpha"))
  en <- filter(summ, grepl("Elastic", model)); bs <- filter(summ, grepl("Baseline", model))
  model_cols <- c("Elastic net, 128 features" = "#1d1d1f", "Baseline, global relative alpha" = "#5f7fb0")
  p6a <- ggplot(roc_df, aes(fpr, tpr, colour = model, group = interaction(model, repeat_id))) +
    geom_abline(linetype = "dashed", colour = scales::alpha(INK, 0.14)) +
    geom_step(alpha = 0.5, linewidth = 0.5) +
    scale_colour_manual(values = model_cols,
      labels = c(sprintf("Baseline, global relative alpha (AUC %.2f)", bs$auc_mean),
                 sprintf("Elastic net, 128 features (AUC %.2f)", en$auc_mean))) +
    coord_equal() +
    labs(x = "False-positive rate", y = "True-positive rate",
         subtitle = "Out-of-fold ROC, one curve per CV repeat\n(10 x repeated 5-fold, nested)") +
    theme_portfolio() + theme(legend.position = "bottom", legend.direction = "vertical",
                              panel.grid.major.x = element_line(colour = scales::alpha(INK, 0.07)))
  p6b <- ggplot(perm, aes(auc)) +
    geom_histogram(binwidth = 0.025, fill = scales::alpha(INK, 0.3), colour = BG, boundary = 0.5) +
    geom_vline(xintercept = en$auc_mean, colour = INK, linewidth = 0.9) +
    annotate("text", x = en$auc_mean, y = Inf, vjust = 1.5, hjust = 1.05, size = 3.4, colour = INK,
             label = sprintf("Observed %.2f\nperm. p = %.3f", en$auc_mean, en$perm_p)) +
    labs(x = "Cross-validated AUC", y = "Permutations",
         subtitle = sprintf("Null distribution, %d label permutations\n(spike at 0.5 = model kept no features)", nrow(perm))) +
    theme_portfolio()
  p6 <- p6a + p6b + plot_layout(widths = c(1, 1.1)) +
    plot_annotation(title = "Band power separates the groups above chance, but far from diagnostically",
                    theme = theme_portfolio())
  save_fig(p6, "fig6_classifier.png", width = 11, height = 5.6)
}
message("Figures written: ", paste(list.files("figures"), collapse = ", "))
