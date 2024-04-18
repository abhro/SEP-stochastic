#!/bin/bash

##SBATCH -J sep
#SBATCH --nodes 1
#####SBATCH --ntasks 5
#SBATCH --cpus-per-task 16
##SBATCH --mem=12000
#SBATCH --mem-per-cpu=10G
#SBATCH --time 03:00:00
#SBATCH --partition=long
#SBATCH --error=logs/%x.%J.err.txt
#SBATCH --output=logs/%x.%J.out.txt

# following trick taken from
# https://arcca.github.io/Introduction-to-Parallel-Programming-using-OpenMP/11-openmp-and-slurm/index.html
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# if environment variables are unset, complain dramatically and die
if [ "x$SHTC_FILE" == x ] || [ "x$MAGGRID_FILE" == x ]; then
  >&2 echo "Please set SHTC_FILE and MAGGRID_FILE"
  exit 1
fi
export SHTC_FILE MAGGRID_FILE

echo "Given parameters:"
echo "SHTC_FILE=$SHTC_FILE"
echo "MAGGRID_FILE=$MAGGRID_FILE"

# this program has a tendency to segfault, so allow
# core dumps and tracebacks in the output
ulimit -c unlimited

echo "Host information:"
echo "hostname: $(hostname)"
echo "hostname -A: $(hostname -A)"
echo "hostname -f: $(hostname -f)"
echo
#echo "Julia is at $(which julia)"
echo "Running maggrid"
./maggrid

if [ $? == 139 ]; then # prev line seg faulted

  echo "Diagnostic info"
  set -v
  echo "Running on node $(hostname -f)"
  cat /proc/sys/kernel/core_pattern


  ls -lR /var/lib/apport

  exit 139
fi
