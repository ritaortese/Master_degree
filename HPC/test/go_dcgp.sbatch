#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --partition=dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH -t 00:10:00
#SBATCH --job-name=stencil
#SBATCH --exclusive

# ===============================

# module load gcc/12.2.0
module load openmpi/4.1.6--gcc--12.2.0

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export TEST_TYPE=${TESTTYPE}
export JOB_NAME=${CSV_NAME}

EXEC=./bin/stencil_parallel

echo "Running with:"
echo "Nodes: $SLURM_JOB_NUM_NODES"
echo "Total tasks: $SLURM_NTASKS"
echo "Tasks per node: $SLURM_TASKS_PER_NODE"
echo "OMP threads: $OMP_NUM_THREADS"
echo "Test type: $TEST_TYPE"
echo "Args: $ARGS"

if [ "$TEST_TYPE" = "omp_strong" ]; then
    ${EXEC} ${ARGS}
else
    mpirun -np ${SLURM_NTASKS} ${EXEC} ${ARGS}
fi