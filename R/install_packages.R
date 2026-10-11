# install_packages.R ---------------------------------------------------------
# Installs the CRAN packages this project uses (skips any already installed).
pkgs <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "scales",
          "patchwork", "lme4", "lmerTest", "glmnet", "pROC", "rmarkdown", "knitr")
missing <- setdiff(pkgs, rownames(installed.packages()))
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, requireNamespace))
message("All packages available.")
