#!/bin/bash

#SBATCH --time=5-00:00:00
#SBATCH --job-name=BLUP
#SBATCH --partition=main
#SBATCH --mem=30G

#SBATCH --array=1-240
#SBATCH --output=/output/dir/BLUP_%a.out
#SBATCH --error=/output/dir/BLUP_%a.err

bash

M=300000

let k=0

for N in 5000 50000; do
for h2 in 0.2 0.5 0.8; do
for poly in 0.001 0.01 0.1 1; do
for n_rep in $(seq 1 10); do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

Rscript /your/blup_analytical_form/pathway/BLUP.R ${N} ${h2} ${poly} ${n_rep} ${M} 

fi
done
done
done
done

