#!/bin/bash

module restore

if [ -n "$WORK_DIR" ]; then
	echo "Changing to $WORK_DIR"
	cd $WORK_DIR
fi

REPO_ROOT=`git rev-parse --show-toplevel`
echo "Script dir = $REPO_ROOT"
source $REPO_ROOT/set_ospray_vars.sh

export LOGFILE="mini-cinema-${SLURM_NNODES}n-${SLURM_JOB_PARTITION}-${SLURM_JOBID}.txt"
export JOB_PREFIX="${SLURM_JOB_PARTITION}-${SLURM_NNODES}n-${SLURM_JOBID}"

ibrun ./mini_cinema \
    -prefix ${JOB_PREFIX} \
    -fif 8 \
    ${JSON_CONFIG} \
    -detailed-stats \
    > ${LOGFILE} 2>&1

