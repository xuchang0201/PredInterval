# PGS-based Phenotype Prediction Interval (PredInterval)

![PredInterval schematic plot](https://github.com/xuchang0201/PredInterval/assets/41645824/f2dc42f0-191c-4566-8ff3-e5a96a2e06ec)

PredInterval is a statistical method that quantifies polygenic score (PGS)-based phenotype prediction uncertainty through the construction of well-calibrated prediction intervals. PredInterval is non-parametric in natural and extracts information based on quantiles of phenotypic residuals through cross-validations, thus achieving well-calibrated coverage of true phenotypic values. In addition, the PredInterval framework is general, takes either individual-level data or summary statistics as input, and can be paired with any PGS method or pre-computed SNP effect sizes obtained from publicly available resources.
    
# How to use PredInterval
There are two versions of PredInterval:
1. **Individual-level version of PredInterval:** construct phenotypic prediction intervals using individual-level genotype and phenotype from the training sample.
2. **Summary statistics version of PredInterval:** construct phenotypic prediction intervals using summary statistics as training data and a small calibration set for phenotypic residual-based calibration.

The PredInterval fitting can be generally divided into three steps:
1. Partitioning training dataset into *k* folds for cross-validation procedure
2. Fitting a PGS method of choice to compute PGSs as inputs for PredInterval
3. Applying PredInterval to construct PGS-based phenotypic prediction intervals

For detailed introduction of model fitting algorithm, please refer to the paper and documentations.  

# Tutorial for Individual-level version of PredInterval
For individual-level version of PredInterval, it requires individual-level genotype and phenotype data of training set and can be fitted based on the following steps:
1. For a pre-specified number of fold *k* (we recommend *k*=5), partition the training set into *k* equal-sized disjoint subsets.
2. For each subset *i* in term, fit a PGS method of choice using the data of the remaining *k*-1 subsets to obtain SNP effect size estimates
3. Apply the SNP effect size estimates from step 2 to compute PGSs for subset *i* and test set using the **score** function in PLINK
4. Repeat step 2 and 3 for *k* times and obtain *k* PGSs for the *k* subsets (e.g., train_PGS_subset_1.profile, ..., train_PGS_subset_k.profile) as well as *k* PGSs for the test set (e.g., test_PGS_subset_1.profile, ..., test_PGS_subset_k.profile)
5. Fit PredInterval to construct phenotypic prediction intervals with the pre-specified confidence level (e.g., 95%)

Example command:
```r
pheno_train=pheno_train.txt
PGS_train_prefix=train_PGS_subset
test_fam=test_set.fam
PGS_test_prefix=test_PGS_subset
cv_fold=5
output=CI_output.txt
conf_level=0.95
Rscript PredInterval.R ${pheno_train} \
${PGS_train_prefix} ${test_fam} \
${PGS_test_prefix} ${cv_fold} ${output} ${conf_level}
```
The inputs and format requirements are:
1. pheno_train: file name for phenotypes of training set. The phenotype file of training set should only include two columns, with first column being sample id and the second column being the phenotypic value.
2. PGS_train_prefix: prefix for the PGSs of *k* subsets of training set from PLINK. Note that the file names of PGSs of *k* training subsets must be named as ${PGS_train_prefix}_subset_1.profile, ..., ${PGS_train_prefix}_subset_k.profile.


# Example
Example codes for fitting PredInterval using toy example to construct 95% confidence interval for PGS-based phenotypic prediction (number of folds=5 for the cross-validation procedure):
1. Fitting individual-level version of PredInterval (toy example data located in the [Individual Level](https://github.com/xuchang0201/PredInterval/tree/main/Toy%20Example/Individual_Level)
```r
workdir=/your/PredInterval/directory
Rscript ${workdir}/PredInterval.R ${workdir}/Individual_Level/pheno_training.txt \
${workdir}/Individual_Level/training_PGS_subset ${workdir}/Individual_Level/test.fam \
${workdir}/Individual_Level/test_PGS_removing_subset 5 ${workdir}/output/CI_ind.txt 0.95
```
2. Fitting summary statistics version of PredInterval (toy example data located in the [Summary Statistics](https://github.com/xuchang0201/PredInterval/tree/main/Toy%20Example/Summary_Statistics)
```r
workdir=/your/PredInterval/directory
Rscript ${workdir}/PredInterval_sumstat.R ${workdir}/Summary_Statistics/pheno_cali_sumstat_demo.txt \
${workdir}/Summary_Statistics/cali_PGS_subset ${workdir}/Summary_Statistics/test_sumstat_demo.fam \
${workdir}/Summary_Statistics/test_sumstat_demo_PGS_removing_subset 5 ${workdir}/output/CI_sumstat.txt 0.95
```



# Summary Statistics from PredInterval Manuscript
To further support reproducibility, we have deposited the summary statistics for the 12 traits in the UKB generated in this manuscript in the google drive link: 
1. [Full Summary statistics](https://drive.google.com/file/d/1Vtb-0IdevPRzhPib-blbjeAgOyUDxj73/view?usp=drive_link): Summary Statistics for the entire training set, used for other benchmark comparison methods
2. [Partitioned Summary statistics](https://drive.google.com/file/d/1IyoqM3ZNTaTMyVCHM1UCjFYCOCtdrw7i/view?usp=drive_link): Summary statistics with each of the five fold removed from the training set, used as input for PredInterval fitting

# Citations

Chang Xu, Santhi K. Ganesh, and Xiang Zhou (2024). Statistical construction of calibrated prediction intervals for polygenic score based phenotype prediction. 

# Questions 
If you have any questions on PredInterval software, please email to Chang Xu (xuchang@umich.edu).
