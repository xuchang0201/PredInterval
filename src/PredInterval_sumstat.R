library(data.table)
args <- commandArgs(TRUE)

# inputs
cali_pheno_name <- args[1] #name of phenotype file for calibration set
pgs_prefix_cali <- args[2] #prefix of PGS output for calibration set 
test_fam_name <- args[3] #name of fam file for test set 
pgs_prefix_test <- args[4] #prefix of PGS output for test set 
n_fold <- as.numeric(args[5]) #number of folds for cross-validation
output_file_name <- args[6] #name of CI output
conf_level <- as.numeric(args[7]) #specify confidence level

# specifying alpha
alpha=1-conf_level

####################### Confidence Interval Construction #######################

### load phenotype file (Y) for calibration set
pheno_cali <- data.frame(fread(cali_pheno_name))
pheno_cali <- pheno_cali[!(is.na(pheno_cali[,2])),] #remove individuals with missing phenotypes in the calibration set 

### load mu(-s(k(i)))(X_i)
mu_calibration <- data.frame(fread(paste0(pgs_prefix_cali, "_", 1, ".profile")))
mu_calibration$subset <- 1
mu_calibration <- mu_calibration[match(pheno_cali$id, mu_calibration$IID), ] #match the order with true phenotype in calibration set

for (k in 2:n_fold){
  pgs_subset_k <- data.frame(fread(paste0(pgs_prefix_cali, "_", k, ".profile")))
  pgs_subset_k <- pgs_subset_k[match(pheno_cali$id, pgs_subset_k$IID), ]  #match the order with true phenotype in calibration set 
  pgs_subset_k$subset <- k
  mu_calibration <- rbind(mu_calibration, pgs_subset_k)
}

### calculate Residual_CV
true_pheno_cali_stacked <- rep(pheno_cali[,2], n_fold)
resid_cv <- abs(true_pheno_cali_stacked-mu_calibration$SCORESUM)
R_i_dataset <- data.frame(mu_calibration$IID, mu_calibration$subset, resid_cv)
colnames(R_i_dataset) <- c("id", "subset", "R_i")

### get result for test set
# test set individuals
test_set <- fread(test_fam_name, header=F)
test_id_list <- test_set$V2

# load PGS for test set with models of different subsets removed
pgs_test_set <- data.frame(fread(paste0(pgs_prefix_test, "_", 1, ".profile")))
pgs_test_set$subset <- 1

for (k in 2:n_fold){
  subsample_kth_test <- data.frame(fread(paste0(pgs_prefix_test, "_", k, ".profile")))
  subsample_kth_test$subset <- k
  pgs_test_set <- rbind(pgs_test_set, subsample_kth_test)
}

### calculate confidence interval for each individual in the test set
lower_bound <- numeric(0)
upper_bound <- numeric(0)

for (j in 1:length(test_id_list)){
  # generate mu(-s(k(i))){X_(n+1)}
  mu_test_set <- rep(pgs_test_set$SCORESUM[which(pgs_test_set$IID==test_id_list[j] & pgs_test_set$subset==1)], length(which(R_i_dataset$subset==1)))
  for (m in 2:n_fold){
    mu_test_set <- c(mu_test_set, rep(pgs_test_set$SCORESUM[which(pgs_test_set$IID==test_id_list[j] & pgs_test_set$subset==m)], length(which(R_i_dataset$subset==m))))
  }
  
  # calculate mu(-s(k(i))){X_(n+1)}-R_i, mu(-s(k(i))){X_(n+1)}+R_i
  mu_minus_R_i <- mu_test_set-R_i_dataset$R_i
  mu_plus_R_i <- mu_test_set+R_i_dataset$R_i
  
  # lower bound
  ordered_mu_minus_R_i <- mu_minus_R_i[order(mu_minus_R_i)]
  lb <- ordered_mu_minus_R_i[ceiling(alpha*(length(ordered_mu_minus_R_i)+1))]
  
  # upper bound
  ordered_mu_plus_R_i <- mu_plus_R_i[order(mu_plus_R_i)]
  ub <- ordered_mu_plus_R_i[ceiling((1-alpha)*(length(ordered_mu_plus_R_i)+1))]
  
  # output
  lower_bound <- append(lower_bound, lb)
  upper_bound <- append(upper_bound, ub)
}

### output confidence interval for all individuals in the test set
output_data <- data.frame(test_id_list, lower_bound, upper_bound)
colnames(output_data)[1] <- "id"
write.table(output_data, output_file_name, row.names=F, col.names=T, quote=F, sep=" ")

