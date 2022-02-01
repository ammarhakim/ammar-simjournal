#!/bin/bash

#SBATCH -A pppl
#SBATCH -n 256
#SBATCH -N 16
#SBATCH -t 24:00:00
#SBATCH -J gkyl
#SBATCH -o gkyl-%j.out
#SBATCH -e gkyl-%j.err
#SBATCH --mail-type=ALL

module load intel
module load intel-mpi

srun --mpi=pmix /home/ammar/gkylsoft/gkyl/bin/gkyl ot1-orsag-tang.lua
