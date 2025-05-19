#!/bin/bash

#SBATCH --time=5-00:00:00
#SBATCH --job-name=LDpred
#SBATCH --partition=main
#SBATCH --mem=30G

#SBATCH --array=1-50
#SBATCH --output=/output/dir/LDpred_%a.out
#SBATCH --error=/output/dir/LDpred_%a.err

bash

rootdir=/your/analysis/rootdir

let k=0

for n_rep in $(seq 1 10); do
for i in $(seq 1 5); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

summstat=/your/sumstat/pathway

nobs=`sed -n "2p" ${summstat} | awk '{print $5}'`
nmis=`sed -n "2p" ${summstat} | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
nsnp=`cat ${summstat} | wc -l`
h2=0.5
ldr=${nsnp}/3000

ldpred_input=/your/ldpred/input/dir

coord=/your/ldpred/coord/dir

ldpred gibbs --cf ${coord}.HDF5 --ldr ${ldr} --ldf /your/LDfile/pathway --out /your/output/pathway --N ${n} --h2 ${h2} --use-gw-h2 --f 0.3 0.1 0.03 0.01 0.003 0.001 0.0003 0.0001

fi
done
done

