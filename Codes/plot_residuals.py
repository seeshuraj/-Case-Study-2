import matplotlib.pyplot as plt
import numpy as np

# Plot for Part 2.2
plt.figure(figsize=(10, 6))
sizes = [8, 16, 32, 64, 128, 256]

for n in sizes:
    residuals = np.loadtxt(f'residuals_n{n}.txt')
    plt.semilogy(residuals, label=f'n={n}')

plt.xlabel('Iteration')
plt.ylabel('Relative Residual (||r_k||/||b||)')
plt.title('GMRES Convergence for Different Problem Sizes')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.savefig('convergence_plot.pdf')
plt.close()

# Plot for Part 2.3 comparing serial vs parallel
plt.figure(figsize=(10, 6))
res_serial = np.loadtxt('residual_serial.txt')
res_parallel = np.loadtxt('residual_parallel.txt')

plt.semilogy(res_serial, label='Serial', linestyle='-')
plt.semilogy(res_parallel, label='Parallel', linestyle='--')
plt.xlabel('Iteration')
plt.ylabel('Relative Residual (||r_k||/||b||)')
plt.title('Serial vs Parallel GMRES Convergence')
plt.legend()
plt.grid(True, which="both", ls="-")
plt.savefig('parallel_comparison.pdf')