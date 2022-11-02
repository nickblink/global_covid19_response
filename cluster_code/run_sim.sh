#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-03:00
#SBATCH -c 10
#SBATCH --mail-type=END
#SBATCH -e errors_outputs/%x_%A_%a.sim.err
#SBATCH -o errors_outputs/%x_%A_%a.sim.out

export R_LIBS_USER="~/R/library"
module load gcc/9.2.0 R/4.2.1 geos/3.10.2 udunits/2.2.28 gdal/3.5.0 jpeg/9b
Rscript R/simulation_main.R $1 $2 $3 $4 ${SLURM_ARRAY_TASK_ID}
