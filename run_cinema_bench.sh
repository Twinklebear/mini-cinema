#!/bin/bash

module restore minicinema-dev

if [ -n "$WORK_DIR" ]; then
	echo "Changing to $WORK_DIR"
	cd $WORK_DIR
fi

REPO_ROOT=`git rev-parse --show-toplevel`
echo "Script dir = $REPO_ROOT"
source $REPO_ROOT/set_ospray_vars.sh

export LOGFILE="mini-cinema-${SLURM_NNODES}n-${SLURM_JOB_PARTITION}-${SLURM_JOBID}.txt"
export JOB_PREFIX="${SLURM_JOB_PARTITION}-${SLURM_NNODES}n-${SLURM_JOBID}"
if [ -n "$LOG_PREFIX" ]; then
    export LOGFILE="${LOG_PREFIX}-${LOGFILE}"
    export JOB_PREFIX="${LOG_PREFIX}-${JOB_PREFIX}"
fi

env

mpirun -n ${SLURM_NNODES} -outfile-pattern mini-cinema-${JOB_PREFIX}-rank-%r.txt \
    ./mini_cinema \
    -prefix $SCRATCH/cinema/${JOB_PREFIX} \
    -fif 1 \
    -detailed-stats \
    ${JSON_CONFIG}

#ibrun ./mini_cinema \
#    -prefix $SCRATCH/cinema/${JOB_PREFIX} \
#    -fif 1 \
#    -detailed-stats \
#    ${JSON_CONFIG} \
#    > ${LOGFILE} 2>&1

