# utils.R --------------------------------------------------------------------
# Shared constants, signal-processing helpers and the plotting theme.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

FS        <- 128                     # sampling rate (Hz)
N_SAMPLES <- 7680                    # 60 s per channel
CHANNELS  <- c("F7", "F3", "F4", "F8", "T3", "C3", "Cz", "C4",
               "T4", "T5", "P3", "Pz", "P4", "T6", "O1", "O2")

BANDS <- tibble::tribble(
  ~band,   ~lo, ~hi,
  "delta",  1,   4,
  "theta",  4,   8,
  "alpha",  8,  13,
  "beta",  13,  30
) |> mutate(band = factor(band, levels = band))
TOTAL_RANGE <- c(1, 30)             # denominator for relative power

# Coarse scalp regions used for the mixed-model analysis
REGIONS <- c(F7 = "Frontal", F3 = "Frontal", F4 = "Frontal", F8 = "Frontal",
             T3 = "Temporal", T4 = "Temporal", T5 = "Temporal", T6 = "Temporal",
             C3 = "Central", Cz = "Central", C4 = "Central",
             P3 = "Parietal", Pz = "Parietal", P4 = "Parietal",
             O1 = "Occipital", O2 = "Occipital")
REGION_LEVELS <- c("Frontal", "Temporal", "Central", "Parietal", "Occipital")

# Approximate 2-D 10-20 positions (nose up), for schematic scalp maps
ELEC_XY <- tibble::tribble(
  ~channel,  ~x,     ~y,
  "F7",    -0.81,  0.59,
  "F3",    -0.40,  0.55,
  "F4",     0.40,  0.55,
  "F8",     0.81,  0.59,
  "T3",    -1.00,  0.00,
  "C3",    -0.50,  0.00,
  "Cz",     0.00,  0.00,
  "C4",     0.50,  0.00,
  "T4",     1.00,  0.00,
  "T5",    -0.81, -0.59,
  "P3",    -0.40, -0.55,
  "Pz",     0.00, -0.50,
  "P4",     0.40, -0.55,
  "T6",     0.81, -0.59,
  "O1",    -0.31, -0.95,
  "O2",     0.31, -0.95
)

# Read one .eea file: a single column of 16 x 7680 values, channel-major.
read_eea <- function(path) {
  v <- scan(path, what = numeric(), quiet = TRUE)
  stopifnot(length(v) == N_SAMPLES * length(CHANNELS))
  m <- matrix(v, nrow = N_SAMPLES, ncol = length(CHANNELS))
  colnames(m) <- CHANNELS
  m
}

# Welch power spectral density for one signal.
# Hann window, linear detrend per segment, one-sided density (units^2 / Hz).
# `keep` optionally gives the start indices of segments to use (artefact rejection).
welch_psd <- function(x, fs = FS, win = 4 * fs, overlap = 0.5, keep = NULL) {
  step   <- floor(win * (1 - overlap))
  starts <- seq(1, length(x) - win + 1, by = step)
  if (!is.null(keep)) starts <- intersect(starts, keep)
  w  <- 0.5 - 0.5 * cos(2 * pi * (0:(win - 1)) / (win - 1))  # Hann
  U  <- sum(w^2)
  tt <- seq_len(win)
  P  <- vapply(starts, function(s) {
    seg <- x[s:(s + win - 1)]
    seg <- resid(lm.fit(cbind(1, tt), seg))                    # linear detrend
    X   <- fft(seg * w)
    Mod(X[1:(win / 2 + 1)])^2 / (fs * U)
  }, numeric(win / 2 + 1))
  p <- rowMeans(P)
  p[2:(win / 2)] <- 2 * p[2:(win / 2)]                         # one-sided
  list(freq = (0:(win / 2)) * fs / win, psd = p, n_segments = length(starts))
}

# Integrate a PSD over [lo, hi) with the rectangle rule.
band_integral <- function(freq, psd, lo, hi) {
  df <- freq[2] - freq[1]
  sum(psd[freq >= lo & freq < hi]) * df
}

# Hedges' g (group2 - group1) with an approximate 95% CI.
hedges_g <- function(x1, x2) {
  n1 <- length(x1); n2 <- length(x2)
  sp <- sqrt(((n1 - 1) * var(x1) + (n2 - 1) * var(x2)) / (n1 + n2 - 2))
  d  <- (mean(x2) - mean(x1)) / sp
  J  <- 1 - 3 / (4 * (n1 + n2) - 9)
  g  <- J * d
  se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
  c(g = g, lo = g - 1.96 * se, hi = g + 1.96 * se)
}

# ---- Plot styling ----------------------------------------------------------
INK      <- "#1d1d1f"
BG       <- "#e8e8ea"
GROUP_COLS <- c(Control = "#607868", Schizophrenia = "#b9844f")
DIVERGING  <- c(low = "#5f7fb0", mid = "#f5f5f7", high = "#b9844f")

theme_portfolio <- function(base_size = 12) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      plot.background   = element_rect(fill = BG, colour = NA),
      panel.background  = element_rect(fill = BG, colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = scales::alpha(INK, 0.07)),
      axis.ticks        = element_blank(),
      axis.text         = element_text(colour = scales::alpha(INK, 0.6)),
      axis.title        = element_text(colour = scales::alpha(INK, 0.7)),
      strip.text        = element_text(colour = INK, face = "bold", size = base_size),
      plot.title        = element_text(colour = INK, face = "bold", size = base_size + 2),
      plot.subtitle     = element_text(colour = scales::alpha(INK, 0.65), size = base_size - 1),
      plot.caption      = element_text(colour = scales::alpha(INK, 0.5), size = base_size - 3, hjust = 0),
      legend.position   = "top",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = scales::alpha(INK, 0.7)),
      plot.margin       = margin(12, 14, 10, 12)
    )
}

save_fig <- function(p, name, width = 8, height = 5) {
  ggsave(file.path("figures", name), p, width = width, height = height,
         dpi = 200, bg = BG)
}
