# PGS-based Phenotype Prediction Interval (PredInterval)

![PredInterval schematic plot](https://github.com/xuchang0201/PredInterval/assets/41645824/f2dc42f0-191c-4566-8ff3-e5a96a2e06ec)

PredInterval is a statistical method that quantifies polygenic score (PGS)-based phenotype prediction uncertainty through the construction of well-calibrated prediction intervals. PredInterval is non-parametric in natural and extracts information based on quantiles of phenotypic residuals through cross-validations, thus achieving well-calibrated coverage of true phenotypic values. In addition, the PredInterval framework is general, takes either individual-level data or summary statistics as input, and can be paired with any PGS method or pre-computed SNP effect sizes obtained from publicly available resources.
    
# How to use PredInterval
There are two versions of PredInterval:
1. **Individual-level version of PredInterval:** construct phenotypic prediction intervals using individual-level genotype and phenotype from the training sample.
2. **Summary statistics version of PredInterval:** construct phenotypic prediction intervals using summary statistics as training data and a small calibration set for phenotypic residual-based calibration.

For detailed introduction of model fitting algorithm, please refer to the paper and documentations.  

# Tutorial for PredInterval
Example codes for the construction of 95% prediction interval for PGS-based phenotype prediction based on 5-fold cross-validations
```r
workdir=/your/data/directory
pheno_train=${workdir}/pheno_train.txt
PGS_train_prefix=${workdir}/train_PGS_subset
test_fam=${workdir}/test_set.fam
PGS_test_prefix=${workdir}/test_PGS_subset
cv_fold=5
output=${workdir}/CI_output.txt
conf_level=0.95
Rscript PredInterval.R ${pheno_train} ${PGS_train_prefix} ${test_fam} ${PGS_test_prefix} ${cv_fold} ${output} ${conf_level} 
```
# Example
Example codes for fitting individual-level version of PredInterval using toy example data:
```r
workdir=/your/PredInterval/directory
Rscript ${workdir}/PredInterval.R ${workdir}/toy_ind/pheno_training.txt \
${workdir}/toy_ind/training_PGS_subset ${workdir}/toy_ind/test.fam \
${workdir}/toy_ind/test_PGS_removing_subset 5 ${workdir}/output/CI_ind.txt 0.95
```

# Summary Statistics from PredInterval Manuscript
To further support reproducibility, we have deposited the summary statistics for the 12 traits in the UKB generated in this manuscript in the google drive link: 
1. [Full Summary statistics](https://drive.google.com/file/d/1Vtb-0IdevPRzhPib-blbjeAgOyUDxj73/view?usp=drive_link): Summary Statistics for the entire training set, used for other benchmark comparison methods
2. [Partitioned Summary statistics](https://drive.google.com/file/d/1IyoqM3ZNTaTMyVCHM1UCjFYCOCtdrw7i/view?usp=drive_link): Summary statistics with each of the five fold removed from the training set, used as input for PredInterval fitting

# Citations

Chang Xu, Santhi K. Ganesh, and Xiang Zhou (2024). Statistical construction of calibrated prediction intervals for polygenic score based phenotype prediction. 

# Questions 
If you have any questions on PredInterval software, please email to Chang Xu (xuchang@umich.edu).
