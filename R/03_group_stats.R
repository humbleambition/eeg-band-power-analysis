# 03_group_stats.R -----------------------------------------------------------
# Group comparisons (Schizophrenia vs Control):
#  1. Mass-univariate: 16 channels x 4 bands, for log absolute and relative
#     power. Welch t-test + Wilcoxon rank-sum (robustness), Hedges' g with 95% CI,
#     Benjamini-Hochberg FDR within each measure (64 tests).
#  2. Region-level linear mixed models per band:
#       log10 power ~ group * region + (1 | subject)
#  3. Individual alpha peak frequency (Welch t-test, Hedges' g).
# Outputs: results/channel_band_stats.csv, results/lmm_anova.csv,
#          results/lmm_group_by_region.csv, results/alpha_peak_stats.csv

source("R/utils.R")
suppressPackageStartupMessages(library(lmerTest))

bp  <- readr::read_csv(file.path("results", "band_power_long.csv"), show_col_types = FALSE) |>
  mutate(group = factor(group, c("Control", "Schizophrenia")),
         band = factor(band, levels(BANDS$band)),
         channel = factor(channel, CHANNELS),
         region = factor(region, REGION_LEVELS))
apf <- readr::read_csv(file.path("results", "alpha_peak.csv"), show_col_types = FALSE) |>
  filter(!is.na(alpha_peak_hz))   # subjects with no identifiable peak are excluded here

# ---- 1. Mass-univariate ----------------------------------------------------
long <- bp |>
  pivot_longer(c(abs_log10, rel), names_to = "measure", values_to = "value")

uni <- long |>
  group_by(measure, band, channel) |>
  summarise({
    x1 <- value[group == "Control"]; x2 <- value[group == "Schizophrenia"]
    tt <- t.test(x2, x1)                                   # Welch
    wt <- wilcox.test(x2, x1, exact = FALSE)
    es <- hedges_g(x1, x2)
    tibble(mean_control = mean(x1), mean_scz = mean(x2),
           t = unname(tt$statistic), df = unname(tt$parameter), p_t = tt$p.value,
           p_wilcox = wt$p.value, g = es[["g"]], g_lo = es[["lo"]], g_hi = es[["hi"]])
  }, .groups = "drop") |>
  group_by(measure) |>
  mutate(p_fdr = p.adjust(p_t, "BH"), p_wilcox_fdr = p.adjust(p_wilcox, "BH")) |>
  ungroup()
readr::write_csv(uni, file.path("results", "channel_band_stats.csv"))

# ---- 2. Linear mixed models by band ----------------------------------------
options(contrasts = c("contr.sum", "contr.poly"))
lmm_anova <- list(); lmm_region <- list()
for (b in levels(bp$band)) {
  d   <- filter(bp, band == b)
  fit <- lmer(abs_log10 ~ group * region + (1 | subject), data = d)
  a   <- anova(fit, type = 3)                               # Satterthwaite df
  lmm_anova[[b]] <- tibble(band = b, term = rownames(a), F = a$`F value`,
                           num_df = a$NumDF, den_df = a$DenDF, p = a$`Pr(>F)`)
  # group difference (SCZ - Control) within each region, from region-averaged data
  lmm_region[[b]] <- d |>
    group_by(subject, group, region) |>
    summarise(v = mean(abs_log10), .groups = "drop") |>
    group_by(region) |>
    summarise({
      x1 <- v[group == "Control"]; x2 <- v[group == "Schizophrenia"]
      es <- hedges_g(x1, x2); tt <- t.test(x2, x1)
      tibble(diff_log10 = mean(x2) - mean(x1), ratio = 10^(mean(x2) - mean(x1)),
             g = es[["g"]], g_lo = es[["lo"]], g_hi = es[["hi"]], p_t = tt$p.value)
    }, .groups = "drop") |>
    mutate(band = b, .before = 1)
}
lmm_anova  <- bind_rows(lmm_anova) |> mutate(p_fdr_group = NA_real_)
grp_rows   <- lmm_anova$term == "group"
lmm_anova$p_fdr_group[grp_rows] <- p.adjust(lmm_anova$p[grp_rows], "BH")
lmm_region <- bind_rows(lmm_region) |> mutate(p_fdr = p.adjust(p_t, "BH"))
readr::write_csv(lmm_anova,  file.path("results", "lmm_anova.csv"))
readr::write_csv(lmm_region, file.path("results", "lmm_group_by_region.csv"))

# ---- 3. Alpha peak frequency -----------------------------------------------
x1 <- apf$alpha_peak_hz[apf$group == "Control"]; x2 <- apf$alpha_peak_hz[apf$group == "Schizophrenia"]
tt <- t.test(x2, x1); es <- hedges_g(x1, x2)
apf_stats <- tibble(n_control = length(x1), n_scz = length(x2), mean_control = mean(x1), sd_control = sd(x1),
                    mean_scz = mean(x2), sd_scz = sd(x2),
                    t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
                    p_wilcox = wilcox.test(x2, x1, exact = FALSE)$p.value,
                    g = es[["g"]], g_lo = es[["lo"]], g_hi = es[["hi"]])
readr::write_csv(apf_stats, file.path("results", "alpha_peak_stats.csv"))

# ---- console summary -------------------------------------------------------
message("Channel x band tests surviving FDR q < .05:")
print(uni |> group_by(measure, band) |>
        summarise(n_sig = sum(p_fdr < .05), median_g = median(g),
                  min_g = min(g), max_g = max(g), .groups = "drop"))
message("LMM group effects:")
print(filter(lmm_anova, term %in% c("group", "group:region")))
message("Alpha peak frequency:"); print(apf_stats)
