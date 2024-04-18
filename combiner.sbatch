#!/bin/bash
##SBATCH -J combiner
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=3G
#SBATCH --time 1:00:00
#SBATCH --partition=eternity
#SBATCH --error=logs/%x.%J.err.txt
#SBATCH --output=logs/%x.%J.out.txt

if [ "x$COMBINER_B1RS_INFILENAME" == x ]; then
  COMBINER_B1RS_INFILENAME=soho/b1rs.dat
fi

if [ "x$COMBINER_MAGGRID_INFILENAME" == x ]; then
  COMBINER_MAGGRID_INFILENAME=soho/maggrid.dat
fi

if [ "x$COMBINER_MAGGRID_OUTFILENAME" == x ]; then
  COMBINER_MAGGRID_OUTFILENAME=soho/maggrid_combined.dat
fi

export COMBINER_B1RS_INFILENAME COMBINER_MAGGRID_INFILENAME COMBINER_MAGGRID_OUTFILENAME

ulimit -c unlimited

./combiner
