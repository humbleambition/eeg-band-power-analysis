# 04_classifier.R ------------------------------------------------------------
# Can resting band power separate the groups out of sample?
#
# Model A: elastic-net logistic regression (glmnet, alpha = 0.5) on 128 features
#          (log absolute + relative power, 16 channels x 4 bands). Alpha peak
#          frequency is left out because it is undefined for 2 subjects.
#          Lambda tuned by inner 5-fold CV inside each outer training fold
#          (nested CV, so the test fold never informs tuning or scaling).
# Model B: one-feature baseline - logistic regression on global relative alpha.
# Evaluation: 10 x repeated stratified 5-fold CV; AUC and balanced accuracy.
# Significance: label-permutation test (B = 200) of the elastic net's AUC.
# Outputs: results/classifier_cv.csv, results/classifier_summary.csv,
#          results/classifier_oof.csv, results/classifier_perm.csv,
#          results/classifier_coefs.csv

source("R/utils.R")
suppressPackageStartupMessages({ library(glmnet); library(pROC) })
set.seed(2026)

N_REPEATS <- 10; K <- 5; N_PERM <- 200

bp  <- readr::read_csv(file.path("results", "band_power_long.csv"), show_col_types = FALSE)

wide <- bp |>
  select(subject, group, channel, band, abs_log10, rel) |>
  pivot_longer(c(abs_log10, rel), names_to = "measure") |>
  mutate(feature = paste(measure, band, channel, sep = "_")) |>
  select(subject, group, feature, value) |>
  pivot_wider(names_from = feature, values_from = value) |>
  arrange(subject)

y <- as.integer(wide$group == "Schizophrenia")
X <- as.matrix(select(wide, -subject, -group))
global_rel_alpha <- rowMeans(X[, grepl("^rel_alpha_", colnames(X))])
message(sprintf("Feature matrix: %d subjects x %d features", nrow(X), ncol(X)))

strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    f[idx] <- sample(rep_len(seq_len(k), length(idx)))
  }
  f
}

fit_predict_enet <- function(Xtr, ytr, Xte) {
  inner <- strat_folds(ytr, 5)
  cv <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.5,
                  foldid = inner, standardize = TRUE)
  list(p = as.numeric(predict(cv, Xte, s = "lambda.1se", type = "response")), fit = cv)
}
fit_predict_base <- function(xtr, ytr, xte) {
  m <- glm(ytr ~ x, family = binomial, data = data.frame(x = xtr, ytr = ytr))
  as.numeric(predict(m, data.frame(x = xte), type = "response"))
}

bal_acc <- function(y, p) {
  yhat <- as.integer(p >= 0.5)
  mean(c(mean(yhat[y == 1] == 1), mean(yhat[y == 0] == 0)))
}
auc_of <- function(y, p) as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<")))

run_cv <- function(y, folds, with_baseline = TRUE, keep_coefs = FALSE) {
  p_en <- p_bs <- numeric(length(y)); coefs <- list()
  for (k in sort(unique(folds))) {
    tr <- folds != k; te <- folds == k
    r <- fit_predict_enet(X[tr, ], y[tr], X[te, , drop = FALSE])
    p_en[te] <- r$p
    if (keep_coefs) coefs[[k]] <- as.matrix(coef(r$fit, s = "lambda.1se"))[-1, 1]
    if (with_baseline) p_bs[te] <- fit_predict_base(global_rel_alpha[tr], y[tr], global_rel_alpha[te])
  }
  list(p_en = p_en, p_bs = p_bs, coefs = coefs)
}

# ---- repeated nested CV ----------------------------------------------------
cv_rows <- list(); oof <- list(); coef_list <- list()
for (r in seq_len(N_REPEATS)) {
  folds <- strat_folds(y, K)
  res <- run_cv(y, folds, keep_coefs = TRUE)
  cv_rows[[r]] <- tibble(
    repeat_id = r,
    model = c("Elastic net (128 features)", "Baseline: global relative alpha"),
    auc = c(auc_of(y, res$p_en), auc_of(y, res$p_bs)),
    balanced_accuracy = c(bal_acc(y, res$p_en), bal_acc(y, res$p_bs)))
  oof[[r]] <- tibble(repeat_id = r, subject = wide$subject, y = y,
                     p_enet = res$p_en, p_baseline = res$p_bs)
  coef_list <- c(coef_list, res$coefs)
}
cv_tbl <- bind_rows(cv_rows)
oof    <- bind_rows(oof)

# coefficient stability: share of the 50 outer fits selecting each feature
coef_mat <- do.call(cbind, coef_list)
coefs <- tibble(feature = rownames(coef_mat),
                selection_freq = rowMeans(coef_mat != 0),
                mean_coef = rowMeans(coef_mat)) |>
  arrange(desc(selection_freq), desc(abs(mean_coef)))

# ---- permutation test ------------------------------------------------------
obs_auc <- mean(cv_tbl$auc[cv_tbl$model == "Elastic net (128 features)"])
perm_auc <- vapply(seq_len(N_PERM), function(b) {
  yp <- sample(y)
  res <- run_cv(yp, strat_folds(yp, K), with_baseline = FALSE)
  auc_of(yp, res$p_en)
}, numeric(1))
p_perm <- (1 + sum(perm_auc >= obs_auc)) / (1 + N_PERM)

summary_tbl <- cv_tbl |>
  group_by(model) |>
  summarise(auc_mean = mean(auc), auc_sd = sd(auc), auc_min = min(auc), auc_max = max(auc),
            bal_acc_mean = mean(balanced_accuracy), bal_acc_sd = sd(balanced_accuracy),
            .groups = "drop") |>
  mutate(perm_p = ifelse(grepl("Elastic", model), p_perm, NA_real_),
         perm_null_mean = ifelse(grepl("Elastic", model), mean(perm_auc), NA_real_),
         perm_null_95 = ifelse(grepl("Elastic", model), unname(quantile(perm_auc, .95)), NA_real_))

readr::write_csv(cv_tbl, file.path("results", "classifier_cv.csv"))
readr::write_csv(summary_tbl, file.path("results", "classifier_summary.csv"))
readr::write_csv(oof, file.path("results", "classifier_oof.csv"))
readr::write_csv(tibble(perm = seq_len(N_PERM), auc = perm_auc), file.path("results", "classifier_perm.csv"))
readr::write_csv(coefs, file.path("results", "classifier_coefs.csv"))
print(summary_tbl, width = Inf)
print(head(coefs, 12))
