# Case Study 2 — Krylov Subspace Methods (MAP55672)

## 📘 Module: MAP55672 (2024–25)
**Topic:** Krylov Subspace Methods  
**Lecturer:** Dr. P. Fritzsch  
**Student:** Seeshuraj Bhoopalan 
**Submission Deadline:** 26th March 2025

---

## 🧠 Overview

This case study involves the step-by-step implementation of the GMRES (Generalized Minimal Residual) algorithm to solve linear systems \( Ax = b \) using Krylov subspace techniques. The tasks include:

1. **Arnoldi Iteration** to construct an orthonormal basis for the Krylov subspace.
2. **Serial GMRES** implementation and performance analysis.
3. **Parallel GMRES** implementation using OpenMP, with a focus on optimizing bottleneck operations.

---

## 📂 Repository Contents

| Codes Folder             | Description                                                        |
|--------------------------|--------------------------------------------------------------------|
| `gmres.cpp`              | Combined implementation: Arnoldi, Serial & Parallel GMRES (OpenMP) |
| `CaseStudy2.pdf`         | Summary report with algorithm descriptions, analysis, and plots    |
| `convergence.pdf` & `parallel_comparison.pdf`       | Semi-log plot of normalized residuals for all matrix sizes         |
| `residuals_n*.txt`       | Output files tracking residual history for each system size        |
| `plot_residuals.py`      | Python script used for generating convergence plot                 |
| `Makefile`               |  Compile automation script for C++ code                            |

---
| Home directory             | Description                                                        |
|--------------------------|--------------------------------------------------------------------|
| `CaseStudy2.pdf`         | Summary report with algorithm descriptions, analysis, and plots    |
| `Codes Folder`       | Combined implementation code of C++ and python with output       |

---
## ⚙️ How to Compile and Run

### 🛠 Requirements
- C++ compiler with OpenMP support (`g++`, `clang++`, etc.)
- Python (for plotting, optional)

### 🧮 Compile the C++ Code and run

```bash
g++ -fopenmp gmres.cpp -o gmres
```
### Run
```bash
./gmres
```
### 🧮 genrate graph python code
```bash
python3 plot_residuals.py
```
