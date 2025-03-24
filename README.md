# Case Study 2 — Krylov Subspace Methods (MAP55672)

## Overview

This repository contains my solution for Case Study 2 of MAP55672 (2024-25), which focuses on Krylov Subspace Methods, specifically the GMRES algorithm.

## Contents

- `gmres.cpp`: Complete C++ implementation including Arnoldi iteration, serial and parallel GMRES using OpenMP.
- `CaseStudy2.pdf`: Summary report with explanations, tables, plots, and output data.
- `convergence.pdf`: Plot showing residual convergence for multiple system sizes.
- `residuals_n*.txt`: Text files containing residuals for each matrix size.

## How to Compile and Run

Make sure you have a C++ compiler with OpenMP support (e.g., `g++`).

### Compile and run:
```bash
g++ -fopenmp parallel_gmres.cpp -o gmres
./gmres
