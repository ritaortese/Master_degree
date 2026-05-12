import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

BASE_PATH = "output/output_parallel"
OUTPUT_DIR = "test_images"

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# =========================
# MPI STRONG SCALING
# =========================
def mpi_strong_plots():
    df = pd.read_csv(os.path.join(BASE_PATH, "mpi_strong.csv"))

    X = df["nodes"].astype(int).values
    times = df["total_time"].values

    order = np.argsort(X)
    X = X[order]
    times = times[order]

    T1 = times[0]
    speedup = T1 / times
    efficiency = speedup / X

    # TIME
    plt.figure()
    plt.plot(X, times, marker='o')
    plt.xlabel("Number of nodes")
    plt.ylabel("Total time (s)")
    plt.title("MPI Strong Scaling - Time")
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "mpi_strong_time.png"))

    # SPEEDUP
    plt.figure()
    plt.plot(X, speedup, marker='o', label="Measured")
    plt.plot(X, X, linestyle='--', label="Ideal")
    plt.xlabel("Number of nodes")
    plt.ylabel("Speedup")
    plt.title("MPI Strong Scaling - Speedup")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "mpi_strong_speedup.png"))

    # EFFICIENCY
    plt.figure()
    plt.plot(X, efficiency, marker='o')
    plt.xlabel("Number of nodes")
    plt.ylabel("Efficiency")
    plt.title("MPI Strong Scaling - Efficiency")
    plt.ylim(0, 1.1)
    plt.yticks(np.arange(0, 1.2, 0.25))
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "mpi_strong_efficiency.png"))


# =========================
# MPI WEAK SCALING
# =========================
def mpi_weak_plots():
    df = pd.read_csv(os.path.join(BASE_PATH, "mpi_weak.csv"))

    X = df["nodes"].astype(int).values
    times = df["total_time"].values
    times_comm = df["comm_time"].values
    times_comp = df["comp_time"].values

    order = np.argsort(X)
    X = X[order]
    times = times[order]
    times_comm = times_comm[order]
    times_comp = times_comp[order]

    T1 = times[0]
    efficiency = T1 / times

    # TIME
    plt.figure()
    plt.plot(X, times, marker='o', label="Total")
    plt.plot(X, times_comm, marker='s', label="Communication")
    plt.plot(X, times_comp, marker='^', label="Computation")
    plt.xlabel("Number of nodes")
    plt.ylabel("Total time (s)")
    plt.title("MPI Weak Scaling - Time")
    plt.ylim(0, 6.25)
    plt.grid()
    plt.legend()
    plt.savefig(os.path.join(OUTPUT_DIR, "mpi_weak_time.png"))

    # EFFICIENCY
    plt.figure()
    plt.plot(X, efficiency, marker='o')
    plt.xlabel("Number of nodes")
    plt.ylabel("Weak Scaling Efficiency")
    plt.title("MPI Weak Scaling - Efficiency")
    plt.ylim(0, 1.1)
    plt.yticks(np.arange(0, 1.2, 0.25))
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "mpi_weak_efficiency.png"))


# =========================
# OPENMP STRONG SCALING
# =========================
def omp_strong_plots():
    df = pd.read_csv(os.path.join(BASE_PATH, "omp_strong.csv"))

    X = df["cpus_per_task"].astype(int).values
    times = df["total_time"].values

    order = np.argsort(X)
    X = X[order]
    times = times[order]

    T1 = times[0]
    speedup = T1 / times
    efficiency = speedup / X

    # TIME
    plt.figure()
    plt.plot(X, times, marker='o')
    plt.xlabel("Number of threads")
    plt.ylabel("Total time (s)")
    plt.title("OpenMP Strong Scaling - Time")
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "omp_strong_time.png"))

    # SPEEDUP
    plt.figure()
    plt.plot(X, speedup, marker='o', label="Measured")
    plt.plot(X, X, linestyle='--', label="Ideal")
    plt.xlabel("Number of threads")
    plt.ylabel("Speedup")
    plt.title("OpenMP Strong Scaling - Speedup")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "omp_strong_speedup.png"))

    # EFFICIENCY
    plt.figure()
    plt.plot(X, efficiency, marker='o')
    plt.xlabel("Number of threads")
    plt.ylabel("Efficiency")
    plt.title("OpenMP Strong Scaling - Efficiency")
    plt.ylim(0, 1.1)
    plt.yticks(np.arange(0, 1.2, 0.4))
    plt.grid()
    plt.savefig(os.path.join(OUTPUT_DIR, "omp_strong_efficiency.png"))


if __name__ == "__main__":

    if os.path.exists(os.path.join(BASE_PATH, "mpi_strong.csv")):
        mpi_strong_plots()

    if os.path.exists(os.path.join(BASE_PATH, "mpi_weak.csv")):
        mpi_weak_plots()

    if os.path.exists(os.path.join(BASE_PATH, "omp_strong.csv")):
        omp_strong_plots()

    print("Plots generated successfully.")