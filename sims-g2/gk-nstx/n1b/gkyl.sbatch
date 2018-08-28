#!/bin/bash

#SBATCH -A pppl
#SBATCH -n 512
#SBATCH -N 32
#SBATCH -t 24:00:00
#SBATCH -J gkyl
#SBATCH -o gkyl-%j.out
#SBATCH -e gkyl-%j.err
#SBATCH --mail-type=ALL

module load intel
module load intel-mpi

srun --mpi=pmix /home/ammar/gkylsoft/gkyl/bin/gkyl n1-Lz4-no-collisions.lua
