#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-00:05
#SBATCH --mail-type=END
#SBATCH -e lm.sim.2.%j.err
#SBATCH -o lm.sim.2.%j.out

mkdir -p output/2
export R_LIBS_USER="~/R/library"
module load gcc/6.2.0 R/3.4.1
Rscript lin-reg-script-2.R $1 $2
