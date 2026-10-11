# 02_features.R --------------------------------------------------------------
# Band-power features per subject x channel x band:
#   abs_log10 : log10 absolute band power (integrated PSD)
#   rel       : band power / total 1-30 Hz power
# plus the individual alpha peak frequency (posterior channels, 7-13 Hz).
# Outputs: results/band_power_long.csv, results/alpha_peak.csv

source("R/utils.R")

psd <- readRDS(file.path("data", "derived", "psd.rds"))

total <- psd |>
  group_by(subject, group, channel) |>
  summarise(total = band_integral(freq, psd, TOTAL_RANGE[1], TOTAL_RANGE[2]),
            .groups = "drop")

bp <- psd |>
  group_by(subject, group, channel) |>
  reframe(band = BANDS$band,
          power = mapply(function(lo, hi) band_integral(freq, psd, lo, hi),
                         BANDS$lo, BANDS$hi)) |>
  left_join(total, by = c("subject", "group", "channel")) |>
  mutate(abs_log10 = log10(power),
         rel       = power / total,
         region    = factor(REGIONS[as.character(channel)], REGION_LEVELS)) |>
  select(subject, group, channel, region, band, power, abs_log10, rel)

# Individual alpha peak frequency: maximum of the mean posterior spectrum
# (O1, O2, P3, P4, Pz) between 7 and 13 Hz, refined with a parabolic fit.
apf <- psd |>
  filter(channel %in% c("O1", "O2", "P3", "P4", "Pz"), freq >= 7, freq <= 13) |>
  group_by(subject, group, freq) |>
  summarise(psd = mean(psd), .groups = "drop") |>
  group_by(subject, group) |>
  summarise(alpha_peak_hz = {
    i <- which.max(psd)
    if (i == 1 || i == length(psd)) {
      NA_real_                     # maximum at the search edge: no identifiable alpha peak
    } else {
      y <- log(psd[(i - 1):(i + 1)]); df <- freq[2] - freq[1]
      freq[i] + df * 0.5 * (y[1] - y[3]) / (y[1] - 2 * y[2] + y[3])
    }
  }, .groups = "drop")

readr::write_csv(bp, file.path("results", "band_power_long.csv"))
readr::write_csv(apf, file.path("results", "alpha_peak.csv"))
message("Subjects without an identifiable alpha peak: ", sum(is.na(apf$alpha_peak_hz)))
message("Band-power rows: ", nrow(bp), "; subjects: ", n_distinct(bp$subject))
