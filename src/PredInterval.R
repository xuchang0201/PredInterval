library(data.table)
args <- commandArgs(TRUE)

# inputs
training_pheno_name <- args[1] #name of phenotype file for training set
pgs_prefix_train    <- args[2] #prefix of PGS output for each subset of training set
test_fam_name       <- args[3] #name of fam file for test set
pgs_prefix_test     <- args[4] #prefix of PGS output for each subset of test set
n_fold              <- as.numeric(args[5]) #number of folds for cross-validation
output_file_name    <- args[6] #name of CI output
conf_level          <- as.numeric(args[7]) #specify confidence level

# specifying alpha
alpha <- 1 - conf_level


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Load optional covariate files for two-stage variance modeling
# Files expected: {pgs_prefix_train}.cvt.txt and {pgs_prefix_test}.cvt.txt
# Format: n x p data.frame with header; first column must be individual ID
# If both files are present, activate heteroscedastic CV+ using K-fold honest procedure:
#   For each fold k, using the (K-1) training folds:
#     Stage 1: fit Y ~ SCORESUM + covariates to get mean model fitted values mu_hat_i
#     Stage 2: fit |resid| ~ SCORESUM + covariates (Gamma GLM) to get fold-k variance model
#   Then on held-out fold k:
#     Compute R_i = |Y_i - mu_hat_i| / sigma_hat_i where mu_hat_i = E[Y|SCORESUM, covariates]
#     (honest: neither mean nor variance model saw fold k individuals during fitting)
#   For test individuals:
#     Predict mu_hat using fold-k mean model, average across K folds
#     Predict sigma_hat using fold-k variance model, average across K folds
# If files are absent, fall back to original CV+ behavior
#
# OPTIONAL: args[8] = path to phenotype file for a subset of test individuals (held-out set H)
#   Format: two columns (id, phenotype), same as training phenotype file
#   When supplied and cvt files are present, activates test held-out mode:
#   For each fold k:
#     - Fit two-stage model on H using fold-k SCORESUM + covariates
#     - Compute nonconformity scores R_i on H using fold-k model
#     - Predict mu_hat and sigma_hat for remaining test individuals using fold-k model
#   Average mu_hat and sigma_hat across K folds for remaining test individuals
#   Construct prediction intervals for remaining test individuals using pooled R_i from H
#   H individuals receive NA in the output (their phenotypes are known)

cvt_train_file     <- paste0(pgs_prefix_train, ".cvt.txt")
cvt_test_file      <- paste0(pgs_prefix_test,  ".cvt.txt")
use_variance_model <- file.exists(cvt_train_file) && file.exists(cvt_test_file)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Parse optional held-out phenotype file argument
use_test_holdout     <- FALSE
holdout_pheno_file <- NULL

if (length(args) >= 8 && nchar(args[8]) > 0 && file.exists(args[8])) {
  holdout_pheno_file <- args[8]
  if (!use_variance_model) {
    warning("Held-out phenotype file supplied but no .cvt.txt files found — ignored.")
  } else {
    use_test_holdout <- TRUE
    cat("Test held-out mode activated.\n")
    cat("  Held-out phenotype file:", holdout_pheno_file, "\n")
  }
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if (use_variance_model) {
  cvt_train <- data.frame(fread(cvt_train_file))
  cvt_test  <- data.frame(fread(cvt_test_file))
  colnames(cvt_train)[1] <- "id"
  colnames(cvt_test)[1]  <- "id"
  
  if (!identical(colnames(cvt_train), colnames(cvt_test))) {
    stop("Column names in train and test covariate files do not match.")
  }
  
  n_cvt        <- ncol(cvt_train) - 1
  cvt_colnames <- colnames(cvt_train)[-1]
  cat("Detected", n_cvt, "covariate(s) in", cvt_train_file, "\n")
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


####################### Confidence Interval Construction #######################

### load phenotype file (Y) for training set
pheno_train <- data.frame(fread(training_pheno_name))
pheno_train <- pheno_train[!(is.na(pheno_train[,2])),] #remove individuals with missing phenotypes in the training set

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Standardize phenotype column names for reliable merging
colnames(pheno_train)[1] <- "id"
colnames(pheno_train)[2] <- "Y"
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### load mu(-s(k(i)))(X_i)
mu_training <- data.frame(fread(paste0(pgs_prefix_train, "_", 1, ".profile")))
mu_training$subset <- 1

for (k in 2:n_fold){
  pgs_subset_k        <- data.frame(fread(paste0(pgs_prefix_train, "_", k, ".profile")))
  pgs_subset_k$subset <- k
  mu_training         <- rbind(mu_training, pgs_subset_k)
}

mu_training <- mu_training[which(mu_training$IID %in% pheno_train$id),] #remove individuals with missing phenotypes in the training set
pheno_train <- pheno_train[match(mu_training$IID, pheno_train$id),]     #match ID between mu_training and phenotype file of training set

### calculate Residual_CV
resid_cv <- abs(pheno_train$Y - mu_training$SCORESUM)

R_i_dataset <- data.frame(pheno_train$id, mu_training$subset, resid_cv)
colnames(R_i_dataset) <- c("id", "subset", "R_i")


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Load test set PGS here (moved up from original position)
# so pgs_test_set is available for test held-out below

### get result for test set
test_set     <- fread(test_fam_name, header = FALSE)
test_id_list <- test_set$V2

pgs_test_set        <- data.frame(fread(paste0(pgs_prefix_test, "_", 1, ".profile")))
pgs_test_set$subset <- 1

for (k in 2:n_fold){
  kth_removed_test        <- data.frame(fread(paste0(pgs_prefix_test, "_", k, ".profile")))
  kth_removed_test$subset <- k
  pgs_test_set            <- rbind(pgs_test_set, kth_removed_test)
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Test held-out mode (for covariate shift setting, where covaraite effect differs between training and test set)
# For each fold k:
#   1. Use fold-k SCORESUM for H individuals + covariates to fit two-stage model on H
#   2. Compute nonconformity scores R_i on H using fold-k model
#   3. Predict mu_hat and sigma_hat for remaining test individuals using fold-k model
# Pool R_i across all K folds (each H individual contributes once per fold ->
# total K * |H| nonconformity scores, one per H-individual per fold)
# Average mu_hat and sigma_hat across K folds for remaining test individuals
# H individuals receive NA in output

if (use_test_holdout) {
  
  # --- Load held-out phenotype file ---
  holdout_pheno <- data.frame(fread(holdout_pheno_file))
  colnames(holdout_pheno)[1] <- "id"
  colnames(holdout_pheno)[2] <- "Y"
  holdout_pheno    <- holdout_pheno[!is.na(holdout_pheno$Y), ]
  holdout_pheno$id <- as.character(holdout_pheno$id)
  
  # Identify H: intersection of supplied IDs with test set and covariate file
  holdout_ids <- intersect(holdout_pheno$id, as.character(test_id_list))
  holdout_ids <- intersect(holdout_ids, as.character(cvt_test$id))
  
  if (length(holdout_ids) == 0) {
    stop("No held-out individuals found in test set with both phenotype and covariate data.")
  }
  
  cat("Test held-out set H: n =", length(holdout_ids), "\n")
  
  # Remaining test individuals (get prediction intervals)
  predict_ids <- setdiff(as.character(test_id_list), holdout_ids)
  cat("Remaining prediction set: n =", length(predict_ids), "\n")
  
  # H covariates and phenotypes
  holdout_cvt <- cvt_test[match(holdout_ids, as.character(cvt_test$id)), ]
  holdout_Y   <- holdout_pheno$Y[match(holdout_ids, holdout_pheno$id)]
  
  # Storage: pool R_i across all K folds for H
  # Each fold contributes |H| nonconformity scores -> total K*|H| scores
  R_i_pool <- numeric(0)
  
  # Storage: mu_hat and sigma_hat for remaining test individuals, one per fold
  mu_hat_predict_mat    <- matrix(NA, nrow = length(predict_ids), ncol = n_fold)
  sigma_hat_predict_mat <- matrix(NA, nrow = length(predict_ids), ncol = n_fold)
  
  mean_formula_h <- as.formula(paste("Y ~ SCORESUM +",         paste(cvt_colnames, collapse = " + ")))
  var_formula_h  <- as.formula(paste("abs_resid ~ SCORESUM +", paste(cvt_colnames, collapse = " + ")))
  
  for (k in 1:n_fold) {
    
    # Build H data frame using fold-k SCORESUM
    holdout_scoresum_k <- pgs_test_set$SCORESUM[
      match(holdout_ids, as.character(pgs_test_set$IID[pgs_test_set$subset == k]))]
    
    holdout_df_k <- data.frame(
      id       = holdout_ids,
      Y        = holdout_Y,
      SCORESUM = holdout_scoresum_k
    )
    holdout_df_k <- merge(holdout_df_k, holdout_cvt, by = "id", all.x = TRUE)
    
    # Stage 1: fit mean model on H using fold-k SCORESUM
    mean_model_h_k        <- lm(mean_formula_h, data = holdout_df_k)
    holdout_df_k$abs_resid <- abs(residuals(mean_model_h_k))
    
    # Stage 2: fit variance model on H using Stage 1 residuals
    var_model_h_k <- glm(var_formula_h, data = holdout_df_k, family = Gamma(link = "log"))
    
    # Compute nonconformity scores on H for fold k
    sigma_hat_h_k <- predict(var_model_h_k, newdata = holdout_df_k, type = "response")
    R_i_k         <- holdout_df_k$abs_resid / sigma_hat_h_k
    
    # Pool R_i across folds
    R_i_pool <- c(R_i_pool, R_i_k)
    
    # Predict mu_hat and sigma_hat for remaining test individuals using fold-k model
    cvt_predict_k          <- cvt_test[match(predict_ids, as.character(cvt_test$id)), ]
    cvt_predict_k$SCORESUM <- pgs_test_set$SCORESUM[
      match(predict_ids, as.character(pgs_test_set$IID[pgs_test_set$subset == k]))]
    
    mu_hat_predict_mat[, k]    <- predict(mean_model_h_k, newdata = cvt_predict_k)
    sigma_hat_predict_mat[, k] <- predict(var_model_h_k,  newdata = cvt_predict_k,
                                          type = "response")
    
    cat("Fold", k, "| H n =", length(holdout_ids),
        "| R_i pool size so far:", length(R_i_pool), "\n")
    
    # Save fold-k models to disk
#    saveRDS(mean_model_h_k, file = paste0(output_file_name, ".holdout_mean_model_fold_", k, ".rds"))
#    saveRDS(var_model_h_k,  file = paste0(output_file_name, ".holdout_var_model_fold_",  k, ".rds"))
  }
  
  # Average mu_hat and sigma_hat across K folds for remaining test individuals
  mu_hat_test    <- setNames(rowMeans(mu_hat_predict_mat),    predict_ids)
  sigma_hat_test <- setNames(rowMeans(sigma_hat_predict_mat), predict_ids)
  
  cat("Total pooled nonconformity scores from H:", length(R_i_pool), "\n")
  
  # Replace R_i_dataset with pooled H residuals
  # subset column set to 1 for all (covariate branch does not use fold structure of R_i)
  R_i_dataset <- data.frame(id = rep(holdout_ids, n_fold),
                            subset = rep(1:n_fold, each = length(holdout_ids)),
                            R_i    = R_i_pool)
  
  # Update test_id_list to prediction set only
  test_id_list <- predict_ids
  
  cat("Test held-out computation complete.\n")
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Standard training-fold-based variance modeling
# Only runs when use_variance_model=TRUE and use_test_holdout=FALSE

if (use_variance_model && !use_test_holdout) {
  
  if (!all(R_i_dataset$id %in% cvt_train$id)) {
    stop("Some training individuals are missing from the train covariate file.")
  }
  
  if (!all(test_id_list %in% cvt_test$id)) {
    stop("Some test individuals are missing from the test covariate file.")
  }
  
  mu_training_uniq <- mu_training[!duplicated(mu_training$IID), ]
  indiv_df <- data.frame(
    id       = mu_training_uniq$IID,
    Y        = pheno_train$Y[match(mu_training_uniq$IID, pheno_train$id)],
    SCORESUM = mu_training_uniq$SCORESUM,
    fold     = mu_training_uniq$subset
  )
  indiv_df <- merge(indiv_df, cvt_train, by = "id", all.x = TRUE)
  
  var_models      <- vector("list", n_fold)
  mean_models     <- vector("list", n_fold)
  sigma_hat_train <- numeric(nrow(indiv_df))
  resid_stage1    <- numeric(nrow(indiv_df))
  
  mean_formula <- as.formula(paste("Y ~ SCORESUM +",         paste(cvt_colnames, collapse = " + ")))
  var_formula  <- as.formula(paste("abs_resid ~ SCORESUM +", paste(cvt_colnames, collapse = " + ")))
  
  for (k in 1:n_fold) {
    
    train_idx <- which(indiv_df$fold != k)
    held_idx  <- which(indiv_df$fold == k)
    train_k   <- indiv_df[train_idx, ]
    held_k    <- indiv_df[held_idx,  ]
    
    mean_model_k      <- lm(mean_formula, data = train_k)
    train_k$abs_resid <- abs(residuals(mean_model_k))
    
    var_model_k      <- glm(var_formula, data = train_k, family = Gamma(link = "log"))
    mean_models[[k]] <- mean_model_k
    var_models[[k]]  <- var_model_k
    
    mu_hat_held            <- predict(mean_model_k, newdata = held_k)
    resid_stage1[held_idx] <- abs(held_k$Y - mu_hat_held)
    
    sigma_hat_held            <- predict(var_model_k, newdata = held_k, type = "response")
    sigma_hat_train[held_idx] <- sigma_hat_held
    
    cat("Fold", k, "data fitted on", length(train_idx), "individuals, with",
        length(held_idx), "individuals withheld.\n")
    
#    saveRDS(mean_model_k, file = paste0(output_file_name, ".mean_model_fold_", k, ".rds"))
#    saveRDS(var_model_k,  file = paste0(output_file_name, ".var_model_fold_",  k, ".rds"))
  }
  
  R_i_dataset$R_i <- resid_stage1[match(R_i_dataset$id, indiv_df$id)] /
    sigma_hat_train[match(R_i_dataset$id, indiv_df$id)]
  
  cvt_test_matched   <- cvt_test[match(test_id_list, cvt_test$id), ]
  mu_hat_test_mat    <- matrix(NA, nrow = length(test_id_list), ncol = n_fold)
  sigma_hat_test_mat <- matrix(NA, nrow = length(test_id_list), ncol = n_fold)
  
  for (k in 1:n_fold) {
    cvt_test_k          <- cvt_test_matched
    cvt_test_k$SCORESUM <- pgs_test_set$SCORESUM[match(test_id_list,
                                                       pgs_test_set$IID[pgs_test_set$subset == k])]
    mu_hat_test_mat[, k]    <- predict(mean_models[[k]], newdata = cvt_test_k)
    sigma_hat_test_mat[, k] <- predict(var_models[[k]],  newdata = cvt_test_k, type = "response")
  }
  
  mu_hat_test    <- setNames(rowMeans(mu_hat_test_mat),    as.character(test_id_list))
  sigma_hat_test <- setNames(rowMeans(sigma_hat_test_mat), as.character(test_id_list))
  
  cat("Variance model: mu_hat and sigma_hat computed for", length(sigma_hat_test),
      "test individuals.\n")
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


### calculate confidence interval for each individual in the test set
lower_bound    <- numeric(length(test_id_list))
upper_bound    <- numeric(length(test_id_list))
predicted_mean <- numeric(length(test_id_list))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Pre-compute all lookup structures outside loop

n_test         <- length(test_id_list)
progress_steps <- floor(seq(1, n_test, length.out = 21))  # report every ~5%

pgs_test_split <- lapply(1:n_fold, function(m) {
  sub_m <- pgs_test_set[pgs_test_set$subset == m, ]
  setNames(sub_m$SCORESUM, as.character(sub_m$IID))
})

fold_sizes <- sapply(1:n_fold, function(m) sum(R_i_dataset$subset == m))

R_i_vec        <- R_i_dataset$R_i
n_train        <- length(R_i_vec)
R_i_order_asc  <- order(R_i_vec)
R_i_sorted_asc <- R_i_vec[R_i_order_asc]
ub_idx         <- ceiling((1 - alpha) * (n_train + 1))

if (use_variance_model) {
  sigma_vec <- setNames(sigma_hat_test[as.character(test_id_list)],
                        as.character(test_id_list))
} else {
  fold_at_ub <- R_i_dataset$subset[R_i_order_asc[ub_idx]]
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cat("Computing prediction intervals for", n_test, "test individuals...\n")

for (j in 1:n_test) {
  
  if (j %in% progress_steps) {
    pct <- round((j / n_test) * 100)
    cat(sprintf("\r  Progress: [%-20s] %3d%%  (%d / %d)",
                paste(rep("=", pct / 5), collapse = ""),
                pct, j, n_test))
    flush.console()
  }
  
  id_j <- as.character(test_id_list[j])
  
  if (use_variance_model) {
    
    mu_j              <- mu_hat_test[id_j]
    sigma_j           <- sigma_hat_test[id_j]
    upper_bound[j]    <- mu_j + sigma_j * R_i_sorted_asc[ub_idx]
    lower_bound[j]    <- mu_j - sigma_j * R_i_sorted_asc[ub_idx]
    predicted_mean[j] <- mu_j
    
  } else {
    
    mu_test_set <- unlist(lapply(1:n_fold, function(m) {
      rep(pgs_test_split[[m]][id_j], fold_sizes[m])
    }))
    
    mu_minus_R_i <- mu_test_set - R_i_vec
    mu_plus_R_i  <- mu_test_set + R_i_vec
    
    ordered_mu_minus_R_i <- mu_minus_R_i[order(mu_minus_R_i)]
    ordered_mu_plus_R_i  <- mu_plus_R_i[order(mu_plus_R_i)]
    
    lower_bound[j]    <- ordered_mu_minus_R_i[ceiling(alpha * (n_train + 1))]
    upper_bound[j]    <- ordered_mu_plus_R_i[ceiling((1 - alpha) * (n_train + 1))]
    predicted_mean[j] <- mean(sapply(1:n_fold, function(m) pgs_test_split[[m]][id_j]))
  }
}
cat("\n  Done.\n")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# MODIFICATION: Output includes all test individuals.
# H individuals (test held-out set) receive NA for all interval columns.
# Remaining test individuals receive their computed intervals as usual.

all_test_ids <- as.character(fread(test_fam_name, header = FALSE)$V2)

output_data <- data.frame(
  id             = all_test_ids,
  lower_bound    = NA_real_,
  upper_bound    = NA_real_,
  predicted_mean = NA_real_
)

# Fill in computed intervals for prediction set individuals
pred_idx <- match(as.character(test_id_list), output_data$id)
output_data$lower_bound[pred_idx]    <- lower_bound
output_data$upper_bound[pred_idx]    <- upper_bound
output_data$predicted_mean[pred_idx] <- predicted_mean
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

### output confidence interval for all individuals in the test set
write.table(output_data, output_file_name,
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")

cat("Output written to", output_file_name, "\n")
cat("  Prediction intervals computed for:", length(test_id_list), "individuals\n")
if (use_test_holdout) {
  cat("  NA (test held-out set):", length(holdout_ids), "individuals\n")
}