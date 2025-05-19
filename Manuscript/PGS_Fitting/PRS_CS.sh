#!/bin/bash

#SBATCH --time=5-00:00:00
#SBATCH --job-name=PRS_CS
#SBATCH --partition=main
#SBATCH --mem=30G

#SBATCH --array=1-50
#SBATCH --output=/output/dir/PRS_CS_%a.out
#SBATCH --error=/output/dir/PRS_CS_%a.err

bash

rootdir=/your/analysis/rootdir
simu_dir=/your/simulation/dir

let k=0

for n_rep in $(seq 1 10); do
for i in $(seq 1 5); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

data=/your/data/pathway
PRSCS_data=/your/prscs_format/pathway

nobs=`sed -n "2p" ${data} | awk '{print $5}'`
nmis=`sed -n "2p" ${data} | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)

awk '{print $2,$6,$7,$9,$11}' ${data} > ${PRSCS_data}
sed -i '1d' ${PRSCS_data}  
sed -i '1i\SNP A1 A2 BETA P' ${PRSCS_data}

export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export OMP_NUM_THREADS=1

for chr in $(seq 1 10)
do
python ${rootdir}/PRScs/PRScs.py --ref_dir=${rootdir}/PRScs/ref_panel/ldblk_1kg_eur --bim_prefix=/target/bim/pathway --sst_file=${PRSCS_data} --n_gwas=${n} --out_dir=/your/output/dir/prefix --chrom=${chr}
done

fi
done
done

