# run_all.R ------------------------------------------------------------------
# Reproduce the whole analysis from scratch (run from the project root):
#   Rscript run_all.R
# Total runtime is roughly 2-3 minutes on a laptop; the download is ~16 MB.

steps <- c("R/00_download_data.R", "R/01_preprocess.R", "R/02_features.R",
           "R/03_group_stats.R", "R/04_classifier.R", "R/05_figures.R")
for (s in steps) {
  message("\n==> ", s)
  source(s, local = new.env(), echo = FALSE)
}

message("\n==> Rendering report")
rmarkdown::render("report/eeg_band_power_report.Rmd",
                  output_format = c("github_document", "html_document"),
                  quiet = TRUE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Done.")
