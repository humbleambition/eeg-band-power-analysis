# 00_download_data.R ---------------------------------------------------------
# Download the MSU adolescent schizophrenia resting-EEG database (~16 MB zipped)
# into data/raw/ and unzip it. Nothing under data/ is committed to git.
#
# Source: Laboratory for Neurophysiology and Neuro-Computer Interfaces,
#         M.V. Lomonosov Moscow State University
#         http://brain.bio.msu.ru/eeg_schizophrenia.htm
# Run from the project root:  Rscript R/00_download_data.R

base_url <- "http://brain.bio.msu.ru/eeg_data/schizophrenia/"
archives <- c(
  norm = "abb4e3e6d50afefd88c1fc50a7179af9",  # 39 healthy adolescents
  sch  = "904096f72c639b14ba3d5e8d8e0bc6b8"   # 45 adolescents with schizophrenia-spectrum disorders
)

raw_dir <- file.path("data", "raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
options(timeout = max(300, getOption("timeout")))

for (grp in names(archives)) {
  zip_path <- file.path(raw_dir, paste0(grp, ".zip"))
  ok <- file.exists(zip_path) && unname(tools::md5sum(zip_path)) == archives[[grp]]
  if (!ok) {
    message("Downloading ", grp, ".zip ...")
    download.file(paste0(base_url, grp, ".zip"), zip_path, mode = "wb", quiet = TRUE)
  }
  md5 <- unname(tools::md5sum(zip_path))
  if (md5 != archives[[grp]]) {
    warning(sprintf("MD5 mismatch for %s.zip (got %s, expected %s). ",
                    grp, md5, archives[[grp]]),
            "The archive may have been updated upstream; results may differ.")
  }
  out_dir <- file.path(raw_dir, grp)
  unlink(out_dir, recursive = TRUE)
  unzip(zip_path, exdir = out_dir)
  n <- length(list.files(out_dir, pattern = "\\.eea$", ignore.case = TRUE, recursive = TRUE))
  message(sprintf("%-4s: %d recordings (md5 %s)", grp, n, md5))
}
