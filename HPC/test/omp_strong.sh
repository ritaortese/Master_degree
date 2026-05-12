#!/bin/bash

SIZES=16384
ITER=500

for T in 1 2 4 8 16 32 64 112
do
    export NODES=1
    export TASKS_PER_NODE=1
    export CPUS_PER_TASK=$T
    export TIME=00:10:00
    export JOBNAME=omp_strong_${T}

    export TESTTYPE=omp_strong
    export CSV_NAME=omp_strong

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