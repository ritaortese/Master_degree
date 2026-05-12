#!/bin/bash

BASE=8192
ITER=500

for P in 1 2 4 8 16
do
    if [ $P -eq 1 ]; then
        export SIZEX=$BASE
        export SIZEY=$BASE
    elif [ $P -eq 2 ]; then
        export SIZEX=$((BASE * 2))
        export SIZEY=$BASE
    elif [ $P -eq 4 ]; then
        export SIZEX=$((BASE * 2))
        export SIZEY=$((BASE * 2))
    elif [ $P -eq 8 ]; then
        export SIZEX=$((BASE * 4))
        export SIZEY=$((BASE * 2))
    elif [ $P -eq 16 ]; then
        export SIZEX=$((BASE * 4))
        export SIZEY=$((BASE * 4))
    fi

    export NODES=$P
    export TASKS_PER_NODE=8
    export CPUS_PER_TASK=14
    export TIME=00:10:00
    export JOBNAME=mpi_weak_${P}

    export TESTTYPE=mpi_weak
    export CSV_NAME=mpi_weak

    export ARGS="-x ${SIZEX} -y ${SIZEY} -n ${ITER}"

    sbatch \
        --nodes=${NODES} \
        --ntasks-per-node=${TASKS_PER_NODE} \
        --cpus-per-task=${CPUS_PER_TASK} \
        --time=${TIME} \
        --job-name=${JOBNAME} \
        --export=ALL,ARGS="${ARGS}" \
        test/go_dcgp.sbatch
done