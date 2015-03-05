#!/bin/bash

#PBS -N s447-5m-lhdi
#-- mail on execution("b"), termination ("e"), or interruption ("a")
#PBS -m aeb
#PBS -M ahakim@pppl.gov
#PBS -l nodes=8:ppn=8
#PBS -l mem=96000mb
#PBS -l walltime=48:00:00
###PBS -q kruskal
#PBS -r n
#PBS -V
#PBS -j oe

NPROCS=`wc -l < $PBS_NODEFILE`
NNODES=`uniq $PBS_NODEFILE | wc -l`

cd $PBS_O_WORKDIR
CMD="/p/gke/software/gkeyll/bin/gkeyll -i s447-5m-lhdi.lua"
mpiexec -np $NPROCS $CMD

exit

