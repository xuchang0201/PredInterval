# PGS-based Phenotype Prediction Interval (PredInterval)

![PredInterval schematic plot](https://github.com/xuchang0201/PredInterval/assets/41645824/f2dc42f0-191c-4566-8ff3-e5a96a2e06ec)

PredInterval is a statistical method that quantifies polygenic score (PGS)-based phenotype prediction uncertainty through the construction of well-calibrated prediction intervals. PredInterval is non-parametric in natural and extracts information based on quantiles of phenotypic residuals through cross-validations, thus achieving well-calibrated coverage of true phenotypic values. In addition, the PredInterval framework is general, takes either individual-level data or summary statistics as input, and can be paired with any PGS method or pre-computed SNP effect sizes obtained from publicly available resources.
    
# How to Use PredInterval
There are two versions of PredInterval:
1. **Non-covariate version of PredInterval:** construct phenotypic prediction intervals using individual-level genotype and phenotype data.
2. **Covariate version of PredInterval:** construct phenotypic prediction intervals using individual-level genotype, phenotype, and covaraites to achieve contextual calibration.
3. **Summary statistics version of PredInterval:** construct phenotypic prediction intervals using summary statistics as training data and a small calibration set for phenotypic residual-based calibration.

The PredInterval fitting can be generally divided into three steps:
1. Partitioning training dataset into *k* folds for cross-validation procedure
2. Fitting a PGS method of choice to compute PGSs as inputs for PredInterval
3. Applying PredInterval to construct PGS-based phenotypic prediction intervals

For detailed introduction of model fitting algorithm, please refer to the paper and documentations.  

# Tutorial for Non-covariate Version of PredInterval
For individual-level non-covariate version of PredInterval, it requires individual-level genotype and phenotype data of training set and can be fitted based on the following steps:
1. For a pre-specified number of fold *k* (we recommend *k*=5), partition the training set into *k* equal-sized disjoint subsets.
2. For each subset *i* in term, fit a PGS method of choice using the data of the remaining *k*-1 subsets to obtain SNP effect size estimates.
3. Apply the SNP effect size estimates from step 2 to compute PGSs for subset *i* and test set using the **score** function in PLINK.
4. Repeat step 2 and 3 for *k* times and obtain *k* PGSs for the *k* subsets (e.g., train_PGS_subset_1.profile, ..., train_PGS_subset_k.profile) as well as *k* PGSs for the test set (e.g., test_PGS_subset_1.profile, ..., test_PGS_subset_k.profile).
5. Fit PredInterval to construct phenotypic prediction intervals with the pre-specified confidence level (e.g., 95%).

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
1. **pheno_train**: file name for phenotypes of training set. The phenotype file of training set should only include two columns, with first column being sample id and the second column being the phenotypic value.
2. **PGS_train_prefix**: prefix for the PGSs of *k* subsets of training set from PLINK. Note that the file names must be named as ${PGS_train_prefix}_subset_1.profile, ${PGS_train_prefix}_subset_2.profile, ..., ${PGS_train_prefix}_subset_k.profile.
3. **test_fam**: fam file of genotype data for the test set. Note that the genotype data of test set should be in PLINK binary format (with bed, bim, and fam files)
4. **PGS_test_prefix**: prefix for the *k* PGSs of test set from PLINK. Similarly, the file names must be named as ${PGS_test_prefix}_subset_1.profile, ${PGS_test_prefix}_subset_2.profile, ..., ${PGS_test_prefix}_subset_k.profile.
5. **cv_fold**: number of folds for cross-validation
6. **output**: file name of output
7. **conf_level**: target confidence level (e.g., 0.95 for 95% confidence level). 

# Tutorial for Covariate Version of PredInterval for Contextual Calibration

## Background: Contextual Calibration and Data Normalization

Contextual calibration refers to the ability of a method to produce calibrated prediction intervals within subgroups defined by covariates (e.g., sex, age, ancestry). Miscalibration within such subgroups can arise when covariate effects are not properly accounted for during data normalization prior.

In many GWAS settings, standard phenotype normalization procedures that account for covariate effects are sufficient to achieve proper calibration. In most cases, removing covariate effects on the phenotypic mean via linear regression is adequate. For example, in typical UK Biobank analyses, effects of sex, age, and age² are adjusted for directly through regression before data analysis.

In less common situations where covariates also influence phenotypic variance (e.g., differing variance between sexes), additional normalization steps may be required. These may include within-stratum quantile normalization or box-cox transformation or the use of heteroscedastic regression models to obtain variance-standardized residuals prior to model fitting.

When such normalization has not been performed, covariates can instead be incorporated directly into the PredInterval fitting step using the covariate-adjusted version described below. In such cases, if covariate effects further differ between the training and test datasets (e.g. if training and test datasets come from two different populations), it is important to use a subset of the test data, rather than the training data, as the hold-out portion during model fitting.

---

## Covariate Version of PredInterval

The covariate version of PredInterval extends the non-covariate version by incorporating covariates within each cross-validation fold. This ensures that covariate effects are accounted for even when the data normalization failed to account for covaraite effects.

Steps 1–4 for fitting the covariate version are identical to the non-covariate version above. The additional requirement is the preparation of a covariate file for both the training and test sets, as described below.

### Step 5: Prepare covariate files

In addition to the inputs required for the non-covariate version, two covariate files must be supplied:

- **Training covariate file**: `${PGS_train_prefix}.cvt.txt`
- **Test covariate file**: `${PGS_test_prefix}.cvt.txt`

These files must be placed in the same folder as the corresponding PGS prefix files. PredInterval will automatically detect them and activate the covariate version. If these files are absent, PredInterval falls back to the standard non-covariate behavior.

#### Format of `.cvt.txt` files

Each file should be a space-delimited plain text file with a header row. The first column must be the sample ID, followed by one column per covariate. The covariate columns must be identical between the training and test files. 

Example `.cvt.txt` file (with sex, age, age squared, genotyping array, and 20 principal components):

```
id sex age age2 array PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20
1000001 1 52 2704 0 0.0312 -0.0021 0.0087 0.0143 -0.0056 0.0201 -0.0034 0.0078 0.0112 -0.0045 0.0067 0.0023 -0.0089 0.0134 0.0056 -0.0012 0.0098 0.0045 -0.0067 0.0123
1000002 0 47 2209 1 0.0289 0.0034 -0.0061 0.0098 0.0073 -0.0145 0.0056 0.0034 -0.0089 0.0112 0.0045 0.0078 0.0023 -0.0056 0.0134 0.0067 -0.0023 0.0089 0.0045 -0.0078
1000003 1 61 3721 0 -0.0103 0.0112 0.0045 -0.0067 0.0156 0.0089 0.0023 -0.0112 0.0067 0.0034 -0.0089 0.0156 0.0045 0.0078 -0.0023 0.0112 0.0034 -0.0056 0.0089 0.0067
```

### Step 6: Fit the covariate version of PredInterval

Once the `.cvt.txt` files are in place, the command is identical to the non-covariate version — no additional arguments are needed:

```bash
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

PredInterval will automatically detect the `.cvt.txt` files and activate the covariate fitting procedure. Progress updates will be printed to the console during model fitting and interval computation.

### Output

The output format is identical to the non-covariate version, with one additional column:

| Column | Description |
|---|---|
| `id` | Sample ID |
| `lower_bound` | Lower bound of the prediction interval |
| `upper_bound` | Upper bound of the prediction interval |
| `predicted_mean` | Predicted phenotypic value averaged across *k* fold-specific mean models |

Additionally, the fitted models for each fold are saved to disk as `.rds` files alongside the output file (e.g., `CI_output.txt.mean_model_fold_1.rds`, `CI_output.txt.var_model_fold_1.rds`) for reuse on separate test datasets without refitting.


# Tutorial for Summary Statistics Version of PredInterval
For summary statistics version of PredInterval, it requires summary statistics of training set and a calibration set as inputs. Specifically, we leverage [PUMAS software](https://github.com/qlu-lab/PUMAS) to mimic the original cross-validation procedure by partitioning the summary statistics of training set into *k* subsampled summary statistics. Details of underlying model and fitting algorithm for PUMAS can be found at PUMAS paper and documentations. 

The summary statistics version of PredInterval can be fitted as follows:
1. Apply PUMAS to obtain *k* subsampled GWAS summary statistics. The example code for partitioning the summary statistics of training set into 5 subsampled summary statistics with 0.8/0.2 ratio is:
```r
Rscript PUMAS.subsampling.R --k 5 --partitions 0.8,0.2 --trait_name ${trait} \
--gwas_path /your/gwas/path/ --ld_path /your/LD/path --output_path /your/output/path/
```
The outputs from the above code are five subsampled summary statistics: ${trait}.gwas.ite1.txt, ${trait}.gwas.ite2.txt, ..., ${trait}.gwas.ite5.txt.

2. For each subsampled summary statistics in turn, fit a PGS method of choice using this data to obtain SNP effect size estimates. 

3. Apply the SNP effect size estimates from step 2 to compute PGSs for calibration and test set using the **score** function in PLINK.
4. Repeat step 2 and 3 for all *k* subsampled summary statistics, and obtain *k* PGSs for the calibration set (e.g., cali_PGS_subset_1.profile, ..., cali_PGS_subset_k.profile) and *k* PGSs for the test set (e.g., test_PGS_subset_1.profile, ..., test_PGS_subset_k.profile).
5. Fit PredInterval to construct phenotypic prediction intervals with the pre-specified confidence level (e.g., 95%).

Example command:
```r
pheno_cali=pheno_cali.txt
PGS_cali_prefix=cali_PGS_subset
test_fam=test_set.fam
PGS_test_prefix=test_PGS_subset
cv_fold=5
output=CI_sumstat.txt
conf_level=0.95
Rscript PredInterval_sumstat.R ${pheno_cali} \
${PGS_cali_prefix} ${test_fam} \
${PGS_test_prefix} 5 ${output} 0.95
```
The inputs and format requirements are:
1. **pheno_cali**: file name for phenotypes of calibration set (two columns only, with first column being sample id and the second column being the phenotypic value).
2. **PGS_cali_prefix**: prefix for the *k* PGSs of calibration set from PLINK. (file names must be named as ${PGS_cali_prefix}_subset_1.profile, ${PGS_cali_prefix}_subset_2.profile, ..., ${PGS_cali_prefix}_subset_k.profile).
3. **test_fam**: fam file of genotype data for the test set (PLINK binary format for genotypes of test set).
4. **PGS_test_prefix**: prefix for the *k* PGSs of test set from PLINK. (file names must be named as ${PGS_test_prefix}_subset_1.profile, ${PGS_test_prefix}_subset_2.profile, ..., ${PGS_test_prefix}_subset_k.profile).
5. **cv_fold**: number of folds for cross-validation
6. **output**: file name of output
7. **conf_level**: target confidence level (e.g., 0.95 for 95% confidence level). 

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
