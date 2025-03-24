### 2.3 Parallel Implementation of GMRES

In this section, a parallel variant of the GMRES algorithm has been implemented using C++ and OpenMP. The main computational bottlenecks in the GMRES algorithm — the dot product, matrix-vector multiplication, and residual norm calculation — were parallelized using `#pragma omp parallel for`.

#### Key Components Parallelized:
- **Inner product (`dot`)**: Parallel reduction across vector elements.
- **Matrix-vector multiplication (`matvec`)**: Parallelized by assigning each row computation to separate threads.
- **Residual norm**: Computed using parallel accumulation of squared differences.

The Arnoldi process and least-squares back-substitution are retained in serial to preserve numerical stability and simplicity. For testing, the algorithm was run with system sizes n = 8, 16, 32, 64, 128, 256 using a fixed iteration count m = n/2, matching the requirement.

The output includes the normalized residual `||r_k|| / ||b||` to evaluate convergence behavior. Residual values are consistent with the serial implementation, verifying correctness.

This implementation demonstrates improved scalability without altering the mathematical integrity of the original GMRES algorithm. A stopping criterion based on a fixed number of iterations was applied for uniformity across tests.

OpenMP's lightweight threading approach proved effective and easy to integrate into the original serial codebase, allowing for enhanced performance on multicore systems.
