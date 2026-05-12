#!/bin/bash

SIZES=16384
ITER=500

for P in 1 2 4 8 16
do
    export NODES=$P
    export TASKS_PER_NODE=8
    export CPUS_PER_TASK=14
    export TIME=00:10:00
    export JOBNAME=mpi_strong_${P}

    export TESTTYPE=mpi_strong
    export CSV_NAME=mpi_strong

    export ARGS="-x ${SIZES} -y ${SIZES} -n ${ITER}"

    sbatch \
          --nodes=${NODES} \
          --ntasks-per-node=${TASKS_PER_NODE} \
          --cpus-per-task=${CPUS_PER_TASK} \
          --time=${TIME} \
          --job-name=${JOBNAME} \
          --export=ALL,ARGS="${ARGS}" \
          test/go_dcgp.sbatch
done