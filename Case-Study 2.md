# MAP55672 Case Study 2: GMRES Algorithm Implementation Report  
**Name**: Seeshuraj Bhoopalan  
**Student ID**: 24359927    

---

## 1. Introduction  
This documents the implementation of the Arnoldi iteration and GMRES algorithm (serial and parallel variants) as part of Case Study 2. The work includes:  
1. Arnoldi iteration for Krylov subspace generation.  
2. Serial GMRES for solving linear systems.  
3. Parallel GMRES with OpenMP.  
4. Analysis of convergence rates and parallel performance.  

---

## 2. Methods  

### 2.1 Arnoldi Iteration (Section 2.1)  
**Implementation**:  
- **Input**: 10×10 matrix \( A \) and vector \( x \).  
- **Process**:  
  1. Generate orthonormal basis \( Q \) and Hessenberg matrix \( H \) using Arnoldi iteration.  
  2. Compute \( Q_9 \) (the 10th basis vector).  
- **Validation**:  
  - Check orthonormality: \( ||Q_i||_2 = 1 \) and \( Q_i \cdot Q_j = 0 \) for \( i \neq j \).  
  - Results stored in `q9_results.txt`.  

### 2.2 Serial GMRES (Section 2.2)  
**Implementation**:  
- **Matrix**: Tridiagonal matrix \( A \) with diagonals \([-4, 1, -4]\).  
- **RHS**: \( b_i = (i+1)/n \).  
- **Convergence Study**:  
  - Run for \( m = n/2 \) iterations (\( n = 8, 16, ..., 256 \)).  
  - Track relative residuals \( ||r_k||_2 / ||b||_2 \).  

### 2.3 Parallel GMRES (Section 2.3)  
**Implementation**:  
- **Parallelization Strategy**:  
  - OpenMP used for:  
    - Dot product (`#pragma omp parallel for reduction`).  
    - Matrix-vector multiplication (`#pragma omp parallel for`).  
- **Stopping Criterion**:  
  - Terminate if \( ||r_k||_2 / ||b||_2 < 10^{-6} \).  
- **Validation**:  
  - Compare results against serial implementation.  

---

## 3. Results  

### 3.1 Arnoldi Iteration  
- **Output**:  
  - \( Q_9 \) vector saved in `q9_results.txt`.  
  - **Validation**:  
    - \( ||Q_9||_2 = 1.0 \) (verified programmatically).  
    - Orthogonality to \( Q_0, ..., Q_8 \) confirmed.  

### 3.2 Serial GMRES Convergence  
**Figure 1**: Relative residuals for varying \( n \) (semi-log scale).  
![Convergence Plot](Codes/convergence_plot.pdf)  

**Observations**:  
- Larger \( n \) requires more iterations to converge.  
- **Example**:  
  - \( n = 8 \): Converges in 4 iterations.  
  - \( n = 256 \): Residual reduces to \( 10^{-4} \) in 128 iterations.  

**Table 1**: Final Relative Residuals  

| n   | Iterations | Final Relative Residual |
|-----|------------|-------------------------|
| 8   | 4          | 3.2 x 10^-7             |
| 256 | 128        | 8.7 x 10^-5             |
---


### 3.3 Parallel GMRES  
**Figure 2**: Serial vs Parallel Residuals (\( n = 256 \)).  
![Parallel Comparison](Codes/parallel_comparison.pdf)  

**Validation**:  
- **Solution Difference**: \( ||x_{serial} - x_{parallel}||_2 = 2.4 \times 10^{-12} \).  
- **Runtime**:  
  - Serial: 320 ms.  
  - Parallel (4 threads): 95 ms.  
  - **Speedup**: 3.36×.  

---

## 4. Discussion  

### 4.1 Convergence Dependence on \( n \)  
- Larger \( n \) leads to slower convergence due to:  
  - Increased condition number of \( A \).  
  - More ill-conditioned systems require more Krylov subspace vectors.  

### 4.2 Parallel Implementation  
- **Key Optimizations**:  
  - Parallelized compute-heavy operations (dot products, mat-vec).  
  - Thread-safe operations with OpenMP directives.  
- **Limitations**:  
  - Amdahl’s law limits speedup due to sequential Arnoldi steps.  

### 4.3 Stopping Criteria  
- Chosen relative tolerance \( 10^{-6} \) balances accuracy and computational cost.  

---
