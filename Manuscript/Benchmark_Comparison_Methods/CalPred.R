devtools::load_all("/your/calpred/pathway")
library(ggplot2); theme_set(theme_classic(base_size=14))
suppressPackageStartupMessages(library(dplyr))
library(patchwork)
library(data.table)
library(dplyr)
args <- commandArgs(TRUE)

#specify inputs
train_fam_file <- args[1] 
test_fam_file <- args[2] 
pheno_all_file <- args[3] 
pgs_train_file <- args[4] 
pgs_test_file <- args[5] 
output_name <- args[6] 

#load covariate file
covar <- fread("/your/covar/file")
covar <- covar[,c("id", "intercept", "sex", "age", "PC1")] 

#load datasets
train_fam <- fread(train_fam_file, header=F)
test_fam <- fread(test_fam_file, header=F)

pheno_all <- fread(pheno_all_file)
pheno_all <- na.omit(pheno_all)
colnames(pheno_all)[2] <- "pheno"
pheno_all <- inner_join(pheno_all, covar, by="id") #merge with covar file


pgs_train <- fread(pgs_train_file)
pgs_test <- fread(pgs_test_file)

#prepare required inputs
#cal_data
pheno_train <- pheno_all[match(train_fam$V2, pheno_all$id), ]
yhat_train <- pgs_train[,c("IID", "SCORESUM")] 
colnames(yhat_train) <- c("id", "yhat")
cal_data <- inner_join(yhat_train, pheno_train, by="id")
cal_data <- subset(cal_data, select=-id)
cal_data <- as.data.frame(cal_data)

#test_data
pheno_test <- pheno_all[match(test_fam$V2, pheno_all$id), ]
yhat_test <- pgs_test[,c("IID", "SCORESUM")]
colnames(yhat_test) <- c("id", "yhat")
test_data <- inner_join(yhat_test, pheno_test, by="id")
test_data <- subset(test_data, select=-id)
test_data <- as.data.frame(test_data)

# build reference to convert quantile <-> value using the calibration dataset
qref <- normalize_reference(cal_data[, 'pheno']) 
cal_data['pheno_q'] <- qref$x2q(cal_data[, "pheno"])
test_data['pheno_q'] <- qref$x2q(test_data[, 'pheno'])

### Use calibration data to train the model
mean_mat_cal <- as.matrix(cal_data[, c("intercept", "yhat", "sex", "age", "PC1")])
sd_mat_cal <- as.matrix(cal_data[, c("intercept", "sex", "age", "PC1"), drop=FALSE])
y_cal <- cal_data[, 'pheno_q']
fit <- train(mean_mat=mean_mat_cal, sd_mat=sd_mat_cal, y=y_cal)

#construct prediction interval
mean_mat_test <- as.matrix(test_data[, c("intercept", "yhat", "sex", "age", "PC1")])
sd_mat_test <- as.matrix(test_data[, c("intercept", "sex", "age", "PC1"), drop=FALSE])
pred <- predict(mean_mat=mean_mat_test, sd_mat=sd_mat_test, mean_coef=fit$mean_coef, sd_coef=fit$sd_coef)
pred <- pred %>%
  rowwise() %>%
  mutate(
    point = qref$q2x(mean), 
    low = qref$q2x(mean - sd * 1.96), 
    high = qref$q2x(mean + sd * 1.96)
  )

#get output
output <- data.frame(test_data, pred)
write.table(output, output_name, row.names=F, col.names=T, quote=F, sep=" ")