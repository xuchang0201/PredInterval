#!/bin/bash

#SBATCH --time=5-00:00:00
#SBATCH --job-name=DBSLMM
#SBATCH --partition=main
#SBATCH --mem=30G

#SBATCH --array=1-600%200
#SBATCH --output=/output/dir/DBSLMM_%a.out
#SBATCH --error=/output/dir/DBSLMM_%a.err

bash

workdir=/your/work/dir

let k=0

for h2 in 0.2 0.5 0.8; do
for poly in 0.001 0.01 0.1 1; do
for n_rep in $(seq 1 10); do
for i in $(seq 1 5); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

data=/your/data/pathway

nobs=`sed -n "2p" ${data} | awk '{print $5}'`
nmis=`sed -n "2p" ${data} | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
nsnp=`cat ${data} | wc -l`

Rscript /your/DBSLMM/pathway/DBSLMM.R --summary ${data} --outPath /your/output/pathway/ --plink /usr/cluster/bin/plink-1.9 --model DBSLMM --dbslmm /your/DBSLMM/pathway/dbslmm --ref /your/refpanel/pathway --n ${n} --nsnp ${nsnp} --block /your/block/pathway/block.txt --h2 ${h2} --thread 5

fi
done
done
done
done
