#!/bin/bash

##SBATCH -J SEProp
#SBATCH --nodes 1
##SBATCH --ntasks 20
#SBATCH --cpus-per-task 20
##SBATCH --mem=12000
#SBATCH --mem=10G
#SBATCH --time 188:00:00
#SBATCH --partition=eternity
#SBATCH --error=logs/%x.%J.err.txt
#SBATCH --output=logs/%x.%J.out.txt

ulimit -s unlimited
export OMP_STACKSIZE=9g
# following trick taken from
# https://arcca.github.io/Introduction-to-Parallel-Programming-using-OpenMP/11-openmp-and-slurm/index.html
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# this program has a tendency to segfault, so allow
# core dumps and tracebacks in the output
ulimit -c unlimited
export PGI_TERM=trace

if [ "x$MAGGRID_FILE" == x ]; then
  MAGGRID_FILE=soho/maggrid.dat
fi
if [ "x$B1RS_FILE" == x ]; then
  B1RS_FILE=soho/b1rs.dat
fi
export MAGGRID_FILE B1RS_FILE
./mapb2s
