#!/bin/bash

#PBS -N s430-is-coal
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m aeb
#PBS -M ahakim@pppl.gov
#PBS -l nodes=16:ppn=8
#PBS -l mem=26000mb
#PBS -l walltime=36:00:00
###PBS -q kruskal
#PBS -r n
#PBS -V
#PBS -j oe

NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

cd $PBS_O_WORKDIR
CMD="/p/gke/software/gkeyll/bin/gkeyll -i s430-is-coal.lua"
mpiexec -np $NPROCS $CMD

exit

