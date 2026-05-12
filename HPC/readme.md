# High Performance Computing project
This project presents the parallelization of a two-dimensional stencil-based heat diffusion problem using MPI and OpenMP. The implementation focuses on distributed domain
decomposition, halo communication using non-blocking MPI operations, and thread-level
parallelism inside each node.
Strong and weak scaling experiments were performed on the Leonardo cluster. The
results show excellent MPI scalability up to 16 nodes, while OpenMP scaling is limited at
high thread counts due to memory bandwidth and synchronization effects.

The complete report is available [here](https://github.com/ritaortese/HPC/blob/f88d1f30934c12ae32b460f41be7227219884008/HPC_report.pdf).