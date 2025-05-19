#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --job-name=CalPred
#SBATCH --partition=main
#SBATCH --mem=20G

#SBATCH --array=1-10
#SBATCH --output=/output/dir/CalPred_%a.out
#SBATCH --error=/output/dir/CalPred_%a.err

bash

workdir=/your/workdir

let k=0

for n_rep in $(seq 1 10); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

Rscript ${workdir}/CalPred.R /your/training/data/fam/file/pathway /your/test/data/fam/file/pathway /your/pheno/file/for/all/individuals/pathway /your/training/PGS/pathway /your/test/PGS/pathway /your/output/pathway  

fi
done
