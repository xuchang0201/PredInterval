library(data.table)
args <- commandArgs(TRUE)

# set parameters
h2 <- as.numeric(args[1])
prevalence <- as.numeric(args[2])
causal_percent <- as.numeric(args[3])
n_rep <- as.numeric(args[4])

# get dimensions
all_samples_geno_bim <- fread("/your/data/directory/all_samples_geno_data.bim", header=F)
all_samples_geno_fam <- fread("/your/data/directory/all_samples_geno_data.fam", header=F)
n_snp=dim(all_samples_geno_bim)[1]
n_sample=dim(all_samples_geno_fam)[1]
simu_dir <- "/your/simulation/directory/" 

#################### simulate SNP effect sizes
if (causal_percent==1){
  #simulate SNP effect sizes
  beta <- rep(0, n_snp)
  for (j in 1:n_snp){
    beta[j] <- rnorm(1, mean=0, sd=sqrt(h2/n_snp))
  }
}else{
  #specify the proportion of causal effects
  n_causal_snps=causal_percent*n_snp
  
  # total index
  total_index <- c(1:n_snp)
  
  #sample the index of causal SNPs that are effective on both traits
  causal_index <- sample(total_index, size=n_causal_snps, replace=F)
  
  ##### simulate betas
  beta <- rep(0, n_snp)
  for (i in 1:length(causal_index)){
    beta[causal_index[i]] <- rnorm(1, mean=0, sd=sqrt(h2/n_causal_snps))
  }
}

####################
maf_data <- fread("/your/data/directory/all_samples_maf.frq")

#calculate non-scaled beta
beta_noscl <- beta/sqrt(2*maf_data$MAF*(1-maf_data$MAF))

#output the simulated beta using PLINK score format
beta_output <- data.frame(maf_data$SNP, maf_data$A1, beta_noscl)
beta_name <- paste0(simu_dir, "simulated_data/", "h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_beta.txt")
write.table(beta_output, beta_name, row.names=F, col.names=F, quote=F, sep=" ")

#calculate X*Beta by PLINK score
score_name <- paste0(simu_dir, "simulated_data/", "h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_X_times_beta")

system(paste0("plink-1.9 --bfile ", "/your/data/directory/all_samples_geno_data",  
              " --score ",  beta_name, " 1 2 3 sum",
              " --out ",    score_name))

x_beta_prod <- fread(paste0(score_name, ".profile"))

#simulate epsilon
epsilon <- rep(0, n_sample)
for (k in 1:n_sample){
  epsilon[k] <- rnorm(1, mean=0, sd=sqrt(1-h2))
}

#calculate the final simulated Y
linear_term <- x_beta_prod$SCORESUM+epsilon
linear_term_scaled <- scale(linear_term) #scale 

#add corresponding intercept to tune case/control ratio
intercept <- log(prevalence/(1-prevalence))
pi <- exp(linear_term_scaled+intercept)/(1+exp(linear_term_scaled+intercept))

y_vector <- numeric(0)

for (i in 1:length(pi)){
  y_i <- rbinom(1,1,pi[i])
  y_vector <- append(y_vector, y_i)
}

#adjust the number of cases to reflect the original case control ratio
if (length(which(y_vector==1)) > prevalence*length(y_vector)){
  original_case_index <- which(y_vector==1)
  case_set_to_0 <- sample(original_case_index, size=length(which(y_vector==1))-prevalence*length(y_vector), replace=F)
  y_vector[case_set_to_0] <- 0
}else{
  original_control_index <- which(y_vector==0)
  control_set_to_1 <- sample(original_control_index, size=prevalence*length(y_vector)-length(which(y_vector==1)), replace=F)
  y_vector[control_set_to_1] <- 1
}

#output simulated phenotypes
pheno_output <- data.frame(x_beta_prod$IID, y_vector)
colnames(pheno_output) <- c("eid", "y")
pheno_output_name <- paste0(simu_dir, "simulated_data/", "h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_binary_pheno_all_samples.txt")
write.table(pheno_output, pheno_output_name, row.names=F, col.names=T, quote=F, sep=" ")

##### prepare GEMMA input phenotype file in the same time
fam <- read.table("/your/data/directory/training_data.fam", header=F, stringsAsFactors=F)

pheno_training <- pheno_output[match(fam$V2, pheno_output$eid),]

pheno_gemma_name <- paste0(simu_dir, "simulated_data/", "h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_binary_pheno_gemma_input.txt")

write.table(pheno_training[,2], pheno_gemma_name, row.names=F, col.names=F, quote=F, sep=" ")

##### run linear model regression by GEMMA
pheno_summstat_output <- paste0("h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_binary_pheno_summstat")

system(paste0("/your/GEMMA/directory/gemma-0.98.1-linux-static",   
              " -bfile ",   "/your/data/directory/training_data",
              " -p ",       pheno_gemma_name, " -lm 1",
              " -outdir ",  paste0(simu_dir, "simulated_data/summstat"),
              " -o ",       pheno_summstat_output))

##### generate summary statistics for remaining training data with the k'th subset removed
for (k in 1:5){
  fam_subset_removed <- fread(paste0("/your/data/directory/", "simu_training_subset", k, "_removed.fam"), header=F)
  pheno_fold <- pheno_output[match(fam_subset_removed$V1, pheno_output$eid),]
  pheno_fold_output_name <- paste0(simu_dir, "simulated_data/", "h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_subset", k, "_removed_binary_pheno.txt")
  write.table(pheno_fold[,2], pheno_fold_output_name, row.names=F, col.names=F, quote=F, sep=" ")
  
  subset_removed_summstat_name <- paste0("h2_", h2, "_prev_", prevalence, "_poly_", causal_percent, "_rep_", n_rep, "_subset", k, "_removed_binary_pheno_summstat")
  
  system(paste0("/your/GEMMA/directory/gemma-0.98.1-linux-static",   
                " -bfile ",   paste0("/your/data/directory/", "simu_training_subset", k, "_removed"),
                " -p ",       pheno_fold_output_name, " -lm 1",
                " -outdir ",  paste0(simu_dir, "simulated_data/summstat"),
                " -o ",       subset_removed_summstat_name))
}



