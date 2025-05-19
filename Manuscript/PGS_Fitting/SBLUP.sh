#!/bin/bash

#SBATCH --time=5-00:00:00
#SBATCH --job-name=SBLUP
#SBATCH --partition=main
#SBATCH --mem=30G

#SBATCH --array=1-50
#SBATCH --output=/output/dir/SBLUP_%a.out
#SBATCH --error=/output/dir/SBLUP_%a.err

bash

rootdir=/your/analysis/root/dir

let k=0

for n_rep in $(seq 1 10); do
for i in $(seq 1 5); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

data=/your/data/pathway
sblup_data=/your/sblup_format/data

nobs=`sed -n "2p" ${data} | awk '{print $5}'`
nmis=`sed -n "2p" ${data} | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)
nsnp=`cat ${data} | wc -l`

awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${data} > ${sblup_data} 
sed -i '1d' ${sblup_data}  
sed -i '1i\SNP A1 A2 freq b se p N' ${sblup_data}
cojo=$(echo "${nsnp}*(1/${h2}-1)" | bc -l)

for j in $(seq 1 10)
do
${rootdir}/gcta_1.93.2beta/gcta64 --bfile /your/ref/panel/pathway --chr ${j} --cojo-file ${sblup_data} --cojo-sblup ${cojo} --cojo-wind 2000 --thread-num 2 --out /your/blup/output/pathway/sblup_beta_chr${j} --diff-freq 0.8
done

fi
done
done

