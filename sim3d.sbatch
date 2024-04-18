#!/bin/bash
##SBATCH -J bk_LC_pfss
#SBATCH --nodes 1
##SBATCH --ntasks 20
#SBATCH --cpus-per-task=20
#SBATCH --mem=6G
#SBATCH --time 288:00:00
#SBATCH --partition=eternity
##SBATCH --partition=long
#SBATCH --error=logs/%x.%J.err.txt
#SBATCH --output=logs/%x.%J.out.txt

# following trick taken from
# https://arcca.github.io/Introduction-to-Parallel-Programming-using-OpenMP/11-openmp-and-slurm/index.html
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# REQUIRES FOLLOWING ENVIRONMENT VARIABLES
#   ROOT_DIR
#   SEEDS_FILE              (allows default)
#   MAGGRID_FILE            (allows default)
#   B1RS_FILE               (allows default)
#   PARAM_OUTDIR_PATH       (allows default)
#   PARAM_NODES             (allows default)
#   CME_DATA_FILE           (allows default)
#   CME_PROP_MODEL_FILE     (allows default)
#   EXTENDED_CME_DATA_FILE  (allows default)

#get the current working directory
echo "This job was submitted from $SLURM_SUBMIT_DIR and I am currently in $PWD"
echo "Starting at date"
echo "Date              = $(date)"
echo "Hostname          = $(hostname -s)"
echo "Working Directory = $PWD"
echo
echo "Number of Nodes Allocated      = $SLURM_JOB_NUM_NODES"
echo "Number of Tasks Allocated      = $SLURM_NTASKS"
echo "Number of Cores/Task Allocated = $SLURM_CPUS_PER_TASK"


# this program has a tendency to segfault, so allow
# core dumps and tracebacks in the output
ulimit -c unlimited

ulimit -s unlimited
#export OMP_STACKSIZE=3g

# prepare ptcl data files

export RUN_DIR="${ROOT_DIR}/run-${SLURM_JOB_ID}"
echo "RUN_DIR = ${RUN_DIR}"

cp -r "$ROOT_DIR/run-0" "$RUN_DIR"

export PTCL_TEMPLATE_FILE="${ROOT_DIR}/loadptcl-template.dat"
export PTCL_FILE="${RUN_DIR}/loadptcl.dat"
#tt1=$(echo "scale = 1; ${RUN_COUNT}*60.0-150.0" | bc -l)
#tt1=$(( ( ${SLURM_JOB_ID}%3 + 1 ) * 60 - 150 ))

# if they didn't specify time to run, throw hissy fit and die
if [ "x$tt1" == x ]; then
  >&2 echo "Please set tt1"
  exit 1
fi
#tt1=1440 # 24 hours = 1440 minutes
# not all nodes have bc
longitude=$( python3 -c "print(round(60 - ${tt1} * 0.009167580))" )
#longitude=$(echo "scale = 1; 60.0 - (${RUN_COUNT}*60.0-150.0) * 0.009167580" | bc -l)
echo "tt1 = $tt1"
echo "longitude = $longitude"
sed "s/!time!/$tt1/; s/!long!/$longitude/" "${PTCL_TEMPLATE_FILE}" > "${PTCL_FILE}"

unset tt1 longitude
# end of preparing ptcl data files


if [ "x$SEEDS_FILE" == x ]; then
  SEEDS_FILE="${RUN_DIR}/seeds.dat"
fi
export SEEDS_FILE

echo "Generating seeds at $SEEDS_FILE"
./seedgen
echo "Finished generating seeds"

if [ "x$MAGGRID_FILE" == x ]; then
  MAGGRID_FILE="${ROOT_DIR}/maggrid.dat"
fi
export MAGGRID_FILE

if [ "x$B1RS_FILE" == x ]; then
  B1RS_FILE="${ROOT_DIR}/b1rs.dat"
fi
export B1RS_FILE

if [ "x$PARAM_OUTDIR_PATH" == x ]; then
  # the trailing slash is important
  PARAM_OUTDIR_PATH="${RUN_DIR}/"
fi

if [ "x$PARAM_NODES" == x ]; then
  PARAM_NODES=20
fi

export PARAM_OUTDIR_PATH PARAM_NODES


# set parameters for cme preparation
if [ "x$CME_DATA_FILE" == x ]; then
  CME_DATA_FILE="${RUN_DIR}/cme.dat"
fi
if [ "x$CME_PROP_MODEL_FILE" == x ]; then
  CME_PROP_MODEL_FILE="${RUN_DIR}/cmepropmodel.dat"
fi
if [ "x$EXTENDED_CME_DATA_FILE" == x ]; then
  EXTENDED_CME_DATA_FILE="${RUN_DIR}/extended_cme.dat"
fi
export CME_DATA_FILE CME_PROP_MODEL_FILE EXTENDED_CME_DATA_FILE

echo "Running sim3d at $(date)"
./sim3d
echo "Finished sim3d at $(date)"
