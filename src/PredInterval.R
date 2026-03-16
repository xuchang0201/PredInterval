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
# If files are absent, fall back to original CV+ behavior with fast pre-sort approximation:
#   sort order of (mu_test_set +/- R_i) approximated by sort order of R_i alone,
#   valid when between-fold SCORESUM variation is small relative to R_i variation

cvt_train_file     <- paste0(pgs_prefix_train, ".cvt.txt")
cvt_test_file      <- paste0(pgs_prefix_test,  ".cvt.txt")
use_variance_model <- file.exists(cvt_train_file) && file.exists(cvt_test_file)

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
# so pgs_test_set is available for sigma_hat_test prediction below

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
# MODIFICATION: For each fold k, use (K-1) training folds to fit a two-stage
# variance model, then standardize R_i for held-out fold k individuals.
#
# Stage 1: fit Y ~ SCORESUM + covariates on (K-1) folds
#          -> compute mu_hat_i = E[Y|SCORESUM, covariates] for held-out fold k
#          -> compute abs residuals |Y_i - mu_hat_i| on (K-1) folds for Stage 2
# Stage 2: fit |resid| ~ SCORESUM + covariates (Gamma GLM, log link) on (K-1) folds
#          -> compute sigma_hat_i for held-out fold k
# Apply:   R_i = |Y_i - mu_hat_i| / sigma_hat_i for held-out fold k individuals
# Test:    predict mu_hat and sigma_hat using fold-k models, average across K folds

if (use_variance_model) {
  
  # Validate: all training individuals must appear in covariate file
  if (!all(R_i_dataset$id %in% cvt_train$id)) {
    stop("Some training individuals are missing from the train covariate file.")
  }
  
  # Validate: all test individuals must appear in test covariate file
  if (!all(test_id_list %in% cvt_test$id)) {
    stop("Some test individuals are missing from the test covariate file.")
  }
  
  # Prepare individual-level training data (one row per individual)
  # Use deduplicated mu_training to avoid duplicate rows per individual
  mu_training_uniq <- mu_training[!duplicated(mu_training$IID), ]
  indiv_df <- data.frame(
    id       = mu_training_uniq$IID,
    Y        = pheno_train$Y[match(mu_training_uniq$IID, pheno_train$id)],
    SCORESUM = mu_training_uniq$SCORESUM,
    fold     = mu_training_uniq$subset
  )
  indiv_df <- merge(indiv_df, cvt_train, by = "id", all.x = TRUE)
  
  # Storage for fold-specific models, per-individual sigma_hat and Stage 1 residuals
  var_models      <- vector("list", n_fold)
  mean_models     <- vector("list", n_fold)
  sigma_hat_train <- numeric(nrow(indiv_df))  # one value per individual
  resid_stage1    <- numeric(nrow(indiv_df))  # |Y_i - mu_hat_i| from Stage 1 mean model
  
  mean_formula <- as.formula(paste("Y ~ SCORESUM +",         paste(cvt_colnames, collapse = " + ")))
  var_formula  <- as.formula(paste("abs_resid ~ SCORESUM +", paste(cvt_colnames, collapse = " + ")))
  
  for (k in 1:n_fold) {
    
    # (K-1) training folds for fold k
    train_idx <- which(indiv_df$fold != k)
    held_idx  <- which(indiv_df$fold == k)
    train_k   <- indiv_df[train_idx, ]
    held_k    <- indiv_df[held_idx,  ]
    
    # Stage 1: fit mean model on (K-1) folds, obtain training absolute residuals
    mean_model_k      <- lm(mean_formula, data = train_k)
    train_k$abs_resid <- abs(residuals(mean_model_k))
    
    # Stage 2: fit variance model on (K-1) folds using Stage 1 absolute residuals
    var_model_k      <- glm(var_formula, data = train_k, family = Gamma(link = "log"))
    mean_models[[k]] <- mean_model_k
    var_models[[k]]  <- var_model_k
    
    # Compute held-out fold k nonconformity scores using Stage 1 fitted values
    # mu_hat_i = E[Y|SCORESUM, covariates]; neither model saw fold k during fitting
    mu_hat_held            <- predict(mean_model_k, newdata = held_k)
    resid_stage1[held_idx] <- abs(held_k$Y - mu_hat_held)
    
    # Compute sigma_hat for held-out fold k using Stage 2 variance model
    sigma_hat_held            <- predict(var_model_k, newdata = held_k, type = "response")
    sigma_hat_train[held_idx] <- sigma_hat_held
    
    cat("Fold", k, "data fitted on", length(train_idx), "individuals, with",
        length(held_idx), "individuals withheld.\n")
    
    # Save fold-k mean and variance models to disk for reuse on separate test data
    saveRDS(mean_model_k, file = paste0(output_file_name, ".mean_model_fold_", k, ".rds"))
    saveRDS(var_model_k,  file = paste0(output_file_name, ".var_model_fold_",  k, ".rds"))
  }
  
  # Replace R_i_dataset$R_i with Stage 1 fitted residuals standardized by sigma_hat
  # (original R_i used raw SCORESUM residuals; now uses mu_hat from mean model)
  R_i_dataset$R_i <- resid_stage1[match(R_i_dataset$id, indiv_df$id)] /
    sigma_hat_train[match(R_i_dataset$id, indiv_df$id)]
  
  # Predict mu_hat and sigma_hat for test individuals across K fold-specific models
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
  
  mu_hat_test    <- rowMeans(mu_hat_test_mat)
  mu_hat_test    <- setNames(mu_hat_test,    as.character(test_id_list))
  sigma_hat_test <- rowMeans(sigma_hat_test_mat)
  sigma_hat_test <- setNames(sigma_hat_test, as.character(test_id_list))
  
  cat("Variance model: mu_hat and sigma_hat computed for", length(sigma_hat_test), "test individuals.\n")
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

# Pre-index pgs_test_set by IID and subset for fast O(1) lookup
pgs_test_split <- lapply(1:n_fold, function(m) {
  sub_m <- pgs_test_set[pgs_test_set$subset == m, ]
  setNames(sub_m$SCORESUM, as.character(sub_m$IID))
})

# Pre-compute fold sizes for R_i_dataset
fold_sizes <- sapply(1:n_fold, function(m) sum(R_i_dataset$subset == m))

# Pre-sort R_i_vec once outside loop.
# Covariate branch:    exact — mu_hat is scalar per individual, sort order of
#                      (mu +/- sigma*R_i) identical to sort order of R_i
# No-covariate branch: approximate — mu_test_set varies by fold but variation
#                      is small relative to R_i variation across n_train individuals;
#                      sort order of (mu_test_set +/- R_i) approximated by R_i order.
#                      fold_at_ub/fold_at_lb give the fold of the quantile-index
#                      individual, used to look up the correct SCORESUM.
# Both bounds use ub_idx because:
#   upper bound: q_{1-alpha}(mu + R_i) -> largest R_i values
#   lower bound: q_{alpha}(mu - R_i)   -> largest R_i values (smallest mu - R_i)
R_i_vec        <- R_i_dataset$R_i
n_train        <- length(R_i_vec)
R_i_order_asc  <- order(R_i_vec)
R_i_sorted_asc <- R_i_vec[R_i_order_asc]
ub_idx         <- ceiling((1 - alpha) * (n_train + 1))

if (use_variance_model) {
  sigma_vec <- setNames(sigma_hat_test[as.character(test_id_list)],
                        as.character(test_id_list))
} else {
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Pre-compute fold membership at quantile index for no-covariate branch:
  # fold_at_ub: fold of the training individual ranked at ub_idx by R_i (for upper bound)
  # fold_at_lb: fold of the training individual ranked at ub_idx by R_i (for lower bound)
  # Both use ub_idx since both bounds correspond to the largest R_i values
  fold_at_ub <- R_i_dataset$subset[R_i_order_asc[ub_idx]]
  fold_at_lb <- R_i_dataset$subset[R_i_order_asc[ub_idx]]  # same index, same fold
  # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    
    # Exact: mu_hat scalar per individual, sigma scalar per individual
    # upper = mu + sigma * q_{1-alpha}(R_i)
    # lower = mu - sigma * q_{1-alpha}(R_i)
    mu_j              <- mean(mu_hat_test_mat[j, ])
    sigma_j           <- sigma_vec[id_j]
    upper_bound[j]    <- mu_j + sigma_j * R_i_sorted_asc[ub_idx]
    lower_bound[j]    <- mu_j - sigma_j * R_i_sorted_asc[ub_idx]
    predicted_mean[j] <- mu_j
    
  } else {
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # use the exact original CV+ for no-covariate branch:
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
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  }
}
cat("\n  Done.\n")

### output confidence interval for all individuals in the test set
output_data <- data.frame(id = test_id_list, lower_bound, upper_bound, predicted_mean)
write.table(output_data, output_file_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = " ")