#!/bin/bash

export NODES=1
export TASKS_PER_NODE=1
export CPUS_PER_TASK=1
export TIME=00:10:00
export JOBNAME=serial_run

export TESTTYPE=serial
export CSV_NAME=serial_results

export ARGS="-x 16384 -y 16384 -n 500"

sbatch \
    --nodes=${NODES} \
    --ntasks-per-node=${TASKS_PER_NODE} \
    --cpus-per-task=${CPUS_PER_TASK} \
    --time=${TIME} \
    --job-name=${JOBNAME} \
    --export=ALL,ARGS="${ARGS}" \
    test/go_dcgp.sbatch