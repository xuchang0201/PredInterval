library(data.table)
library(genio)
library(dplyr)
args <- commandArgs(TRUE)

N <- as.numeric(args[1]) #sample size
h2 <- as.numeric(args[2])
poly <- as.numeric(args[3])
n_rep <- as.numeric(args[4])
M <- as.numeric(args[5]) #number of SNPs

# load genotype matrix for testing set in the simulations
bim <- read_bim("/your/genotype/data/pathway/testing_data.bim")
fam <- read_fam("/your/genotype/data/pathway/testing_data.fam")
test_geno <- read_bed("/your/genotype/data/pathway/testing_data.bed", bim$id, fam$id)
test_geno <- t(test_geno) #N times M matrix

#impute missing SNPs
Z <- matrix(0, ncol=ncol(test_geno), nrow=nrow(test_geno))
for (j in 1:ncol(test_geno)){
  Z[,j] <- ifelse(is.na(test_geno[,j]), mean(test_geno[,j], na.rm=TRUE), test_geno[,j])
}

#standardize genotype matrix
geno_s <- scale(Z, center=T, scale=T)

#compute vectors of X^2
x_square_vec <- numeric(0)

for (k in 1:nrow(geno_s)){
  x_ij_square_sum <- sum((geno_s[k,])^2)
  x_square_vec <- append(x_square_vec, x_ij_square_sum)
}

#output the result
sum_x_sq <- data.frame(eid=row.names(test_geno), sum_x_sq=x_square_vec)

###### BLUP analytical form computation
#load PGS point estimates from SBLUP
pgs <- data.frame(fread("/your/PGS/file/pathway"))
pgs <- pgs[,c("IID", "SCORESUM")]
colnames(pgs) <- c("eid", "point")

#compute CI
ci_result <- inner_join(pgs, sum_x_sq, by="eid")

ci_result$var <- (ci_result$sum_x_sq/N)*(1/(1+M/(h2*N)))
ci_result$low <- ci_result$point-1.96*sqrt(ci_result$var+(1-h2))
ci_result$high <- ci_result$point+1.96*sqrt(ci_result$var+(1-h2))

#merge with true pheno in the test set
pheno_all <- data.frame(fread("/your/pheno/file/for/all/individuals/pathway"))
colnames(pheno_all)[2] <- "pheno"
output <- inner_join(ci_result, pheno_all, by="eid")

#output
write.table(output, "/your/output/pathway", row.names=F, col.names=T, quote=F, sep=" ")