# 01_preprocess.R ------------------------------------------------------------
# Load every recording, run quality checks, reject artefactual 4-s segments,
# and estimate a Welch PSD per subject x channel.
# Outputs:
#   data/derived/psd.rds          long PSD table (subject, group, channel, freq, psd)
#   results/qc_subject.csv        per-subject QC summary
#   results/qc_channel_flags.csv  channels with outlying variance (robust z > 3.5)

source("R/utils.R")
set.seed(1)

raw_dir <- file.path("data", "raw")
dir.create(file.path("data", "derived"), recursive = TRUE, showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

files <- tibble(path = list.files(raw_dir, pattern = "\\.eea$", recursive = TRUE,
                                  full.names = TRUE, ignore.case = TRUE)) |>
  mutate(group   = ifelse(grepl("/norm/", path), "Control", "Schizophrenia"),
         subject = paste0(ifelse(group == "Control", "C_", "S_"),
                          toupper(sub("\\.eea$", "", basename(path), ignore.case = TRUE))))
stopifnot(!anyDuplicated(files$subject))
message("Recordings found: ", nrow(files), " (",
        paste(names(table(files$group)), table(files$group), collapse = ", "), ")")

WIN <- 4 * FS; STEP <- WIN / 2
seg_starts <- seq(1, N_SAMPLES - WIN + 1, by = STEP)          # 29 overlapping 4-s segments

process_one <- function(path, subject, group) {
  m <- read_eea(path)
  m <- sweep(m, 2, colMeans(m))                                 # remove DC offset

  # --- channel-level checks ---------------------------------------------------
  ch_sd   <- apply(m, 2, sd)
  flat    <- ch_sd < 1e-6
  n_na    <- sum(is.na(m))

  # --- segment-level artefact rejection --------------------------------------
  # Peak-to-peak amplitude per segment and channel; a segment is rejected when
  # any channel exceeds that channel's median + 5 * MAD (within subject).
  ptp <- sapply(seg_starts, function(s) apply(m[s:(s + WIN - 1), ], 2, function(z) diff(range(z))))
  thr <- apply(ptp, 1, function(z) median(z) + 5 * mad(z))
  bad_seg <- apply(ptp > thr, 2, any)
  keep <- seg_starts[!bad_seg]

  # --- spectra ---------------------------------------------------------------
  spec <- lapply(CHANNELS, function(ch) {
    w <- welch_psd(m[, ch], keep = keep)
    tibble(channel = ch, freq = w$freq, psd = w$psd)
  }) |> bind_rows()

  # line-noise check at 50 Hz (ratio of 49-51 Hz to flanking 44-48 / 52-56 Hz)
  ln <- spec |>
    group_by(channel) |>
    summarise(r = mean(psd[freq >= 49 & freq <= 51]) /
                  mean(psd[(freq >= 44 & freq <= 48) | (freq >= 52 & freq <= 56)]),
              .groups = "drop")

  list(
    psd = spec |> mutate(subject = subject, group = group),
    qc  = tibble(subject = subject, group = group, n_na = n_na,
                 n_flat_channels = sum(flat),
                 segments_total = length(seg_starts),
                 segments_rejected = sum(bad_seg),
                 median_channel_sd = median(ch_sd),
                 line_noise_ratio_50hz = median(ln$r)),
    ch_logvar = tibble(subject = subject, group = group, channel = CHANNELS,
                       logvar = log(ch_sd^2))
  )
}

out <- purrr::pmap(files[, c("path", "subject", "group")], process_one)

psd <- bind_rows(lapply(out, `[[`, "psd")) |>
  mutate(group = factor(group, c("Control", "Schizophrenia")),
         channel = factor(channel, CHANNELS))
qc  <- bind_rows(lapply(out, `[[`, "qc"))

# Channels with outlying overall variance relative to the same channel in all
# other subjects (robust z on log variance). Reported; not removed, because the
# database was supplied as artefact-free and log band power is robust to gain.
ch_flags <- bind_rows(lapply(out, `[[`, "ch_logvar")) |>
  group_by(channel) |>
  mutate(robust_z = (logvar - median(logvar)) / mad(logvar)) |>
  ungroup() |>
  filter(abs(robust_z) > 3.5) |>
  arrange(desc(abs(robust_z)))

# Exclusion rule (pre-specified): drop a subject if >25% of segments rejected
# or any channel is flat / missing.
qc <- qc |>
  mutate(excluded = segments_rejected / segments_total > 0.25 |
                    n_flat_channels > 0 | n_na > 0)

saveRDS(psd |> filter(!subject %in% qc$subject[qc$excluded]),
        file.path("data", "derived", "psd.rds"))
readr::write_csv(qc, file.path("results", "qc_subject.csv"))
readr::write_csv(ch_flags, file.path("results", "qc_channel_flags.csv"))

message(sprintf("Segments rejected: median %d / %d per subject (max %d).",
                median(qc$segments_rejected), length(seg_starts), max(qc$segments_rejected)))
message("Subjects excluded: ", sum(qc$excluded))
message("Channel variance outliers (|robust z| > 3.5): ", nrow(ch_flags))
message(sprintf("Median 50-Hz line-noise ratio: %.2f", median(qc$line_noise_ratio_50hz)))
