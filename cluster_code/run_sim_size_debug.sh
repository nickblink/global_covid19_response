#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-10:00
#SBATCH -c 10
#SBATCH --mem=40G
#SBATCH --mail-type=END
#SBATCH -x compute-f-17-[09-25]
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

module load gcc/9.2.0 R/4.2.1 geos/3.10.2 udunits/2.2.28 gdal/3.5.0 jpeg/9b cmake/3.22.2
Rscript R/simulation_main_size_debug.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${SLURM_ARRAY_TASK_ID}
