#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <omp.h>

// TYPE DEFINITIONS
typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;

// HELPER FUNCTIONS

// Compute dot product of two vectors
double dot_product(const Vector& a, const Vector& b) {
    double result = 0.0;
    #pragma omp parallel for reduction(+:result)
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}

// Add two vectors
Vector vector_add(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Subtract two vectors
Vector vector_subtract(const Vector& a, const Vector& b) {
    Vector result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// Multiply vector by scalar
Vector vector_scale(const Vector& v, double scalar) {
    Vector result(v.size());
    for (size_t i = 0; i < v.size(); ++i) {
        result[i] = v[i] * scalar;
    }
    return result;
}

// Compute vector norm
double vector_norm(const Vector& v) {
    return std::sqrt(dot_product(v, v));
}

// Matrix-vector multiplication
Vector matrix_vector_mult(const Matrix& A, const Vector& x) {
    Vector result(A.size(), 0.0);
    #pragma omp parallel for
    for (size_t i = 0; i < A.size(); ++i) {
        double temp = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j) {
            temp += A[i][j] * x[j];
        }
        result[i] = temp;
    }
    return result;
}

// ARNOLDI ITERATION
void arnoldi_iteration(const Matrix& A, const Vector& u, int m, Matrix& Q, Matrix& H) {
    int n = A.size();
    Q.resize(m + 1, Vector(n, 0.0));
    H.resize(m + 1, Vector(m, 0.0));

    // Normalize initial vector
    Vector q0 = u;
    double beta = vector_norm(q0);
    if (beta < 1e-14) return;
    q0 = vector_scale(q0, 1.0 / beta);
    Q[0] = q0;

    for (int j = 0; j < m; ++j) {
        Vector v = matrix_vector_mult(A, Q[j]);

        // Orthogonalize against previous vectors
        for (int i = 0; i <= j; ++i) {
            H[i][j] = dot_product(Q[i], v);
            v = vector_subtract(v, vector_scale(Q[i], H[i][j]));
        }

        H[j+1][j] = vector_norm(v);
        
        if (H[j+1][j] > 1e-14) {
            Q[j+1] = vector_scale(v, 1.0 / H[j+1][j]);
        } else {
            break;  // Early termination
        }
    }
}

// GMRES SOLVER 
Vector gmres_solver(const Matrix& A, const Vector& b, int m, Vector& residuals, double tol = 1e-6) {
    int n = A.size();
    Vector x(n, 0.0);  // Initial guess
    
    // Compute initial residual
    Vector r = vector_subtract(b, matrix_vector_mult(A, x));
    double beta = vector_norm(r);
    residuals.clear();
    residuals.push_back(beta);
    
    if (beta < 1e-14) return x;
    double rel_tol = beta * tol;

    Matrix Q(m + 1, Vector(n, 0.0));
    Matrix H(m + 1, Vector(m, 0.0));
    Vector c(m, 0.0), s(m, 0.0), g(m + 1, 0.0);
    g[0] = beta;

    Q[0] = vector_scale(r, 1.0 / beta);

    for (int j = 0; j < m; ++j) {
        // Arnoldi step
        Vector v = matrix_vector_mult(A, Q[j]);
        for (int i = 0; i <= j; ++i) {
            H[i][j] = dot_product(Q[i], v);
            v = vector_subtract(v, vector_scale(Q[i], H[i][j]));
        }
        H[j+1][j] = vector_norm(v);
        
        if (H[j+1][j] > 1e-14) {
            Q[j+1] = vector_scale(v, 1.0 / H[j+1][j]);
        } else {
            break;
        }

        // Apply previous Givens rotations
        for (int i = 0; i < j; ++i) {
            double temp = c[i] * H[i][j] + s[i] * H[i+1][j];
            H[i+1][j] = -s[i] * H[i][j] + c[i] * H[i+1][j];
            H[i][j] = temp;
        }

        // Compute new Givens rotation
        double nu = std::sqrt(H[j][j] * H[j][j] + H[j+1][j] * H[j+1][j]);
        if (nu == 0) {
            c[j] = 1.0;
            s[j] = 0.0;
        } else {
            c[j] = H[j][j] / nu;
            s[j] = H[j+1][j] / nu;
        }

        H[j][j] = nu;
        H[j+1][j] = 0.0;

        double g_j = c[j] * g[j] + s[j] * g[j+1];
        g[j+1] = -s[j] * g[j] + c[j] * g[j+1];
        g[j] = g_j;

        residuals.push_back(std::abs(g[j+1]));

        if (std::abs(g[j+1]) < rel_tol) {
            break;
        }
    }

    // Solve upper triangular system
    int steps = std::min(m, static_cast<int>(residuals.size()) - 1);
    Vector y(steps, 0.0);
    for (int i = steps - 1; i >= 0; --i) {
        y[i] = g[i];
        for (int k = i + 1; k < steps; ++k) {
            y[i] -= H[i][k] * y[k];
        }
        y[i] /= H[i][i];
    }

    // Form solution
    for (int i = 0; i < steps; ++i) {
        x = vector_add(x, vector_scale(Q[i], y[i]));
    }

    return x;
}

// MATRIX GENERATION 
Matrix create_tridiagonal_matrix(int n) {
    Matrix A(n, Vector(n, 0.0));
    for (int i = 0; i < n; ++i) {
        A[i][i] = -4.0;
        if (i > 0) A[i][i-1] = 1.0;
        if (i < n-1) A[i][i+1] = 1.0;
    }
    return A;
}

Vector create_rhs_vector(int n) {
    Vector b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (i + 1.0) / n;
    }
    return b;
}

// TEST CASES
void part2_1() {
    Matrix A = {
        {3,8,7,3,3,7,2,3,4,8},
        {5,4,1,6,9,8,3,7,1,9},
        {3,6,9,4,8,6,5,6,6,6},
        {5,3,4,7,4,9,2,3,5,1},
        {4,4,2,1,7,4,2,2,4,5},
        {4,2,8,6,6,5,2,1,1,2},
        {2,8,9,5,2,9,4,7,3,3},
        {9,3,2,2,7,3,4,8,7,7},
        {9,1,9,3,3,1,2,7,7,1},
        {9,3,2,2,6,4,4,7,3,5}
    };
    
    Vector x = {
        0.757516242460009,
        2.73405796361329,
        -0.555605907443403,
        1.144284746786790,
        0.645280108318073,
        -0.085488474462339,
        -0.623679022063185,
        -0.465240896342741,
        2.382909057772335,
        -0.120465395885881
    };

    Matrix Q;
    Matrix H;
    
    arnoldi_iteration(A, x, 9, Q, H);
    
    std::ofstream out_file("q9_results.txt");
    out_file << "Q_9 vector:\n";
    for (double val : Q[9]) {
        out_file << std::setprecision(15) << val << "\n";
    }
    out_file.close();
}

void part2_2() {
    std::vector<int> problem_sizes = {8, 16, 32, 64, 128, 256};
    
    for (int n : problem_sizes) {
        Matrix A = create_tridiagonal_matrix(n);
        Vector b = create_rhs_vector(n);
        Vector residuals;
        
        Vector x = gmres_solver(A, b, n/2, residuals);
        
        // Save residuals
        double b_norm = vector_norm(b);
        std::ofstream res_file("residuals_n" + std::to_string(n) + ".txt");
        for (double r : residuals) {
            res_file << r / b_norm << "\n";
        }
        res_file.close();
    }
}

void part2_3() {
    int n = 256;
    Matrix A = create_tridiagonal_matrix(n);
    Vector b = create_rhs_vector(n);
    
    // Serial version
    Vector residuals_serial;
    Vector x_serial = gmres_solver(A, b, n/2, residuals_serial);
    
    // Parallel version
    Vector residuals_parallel;
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "Running with " << omp_get_num_threads() << " threads\n";
    }
    Vector x_parallel = gmres_solver(A, b, n/2, residuals_parallel);
    
    // Verify results
    double difference = vector_norm(vector_subtract(x_serial, x_parallel));
    std::cout << "Solution difference: " << difference << "\n";
    
    // Save residuals
    double b_norm = vector_norm(b);
    std::ofstream ser_file("residual_serial.txt");
    std::ofstream par_file("residual_parallel.txt");
    
    for (double r : residuals_serial) ser_file << r / b_norm << "\n";
    for (double r : residuals_parallel) par_file << r / b_norm << "\n";
    
    ser_file.close();
    par_file.close();
}

// MAIN PROGRAM
int main() {
    std::cout << "Running Part 2.1 (Arnoldi Iteration)...\n";
    part2_1();
    
    std::cout << "\nRunning Part 2.2 (Serial GMRES)...\n";
    part2_2();
    
    std::cout << "\nRunning Part 2.3 (Parallel GMRES)...\n";
    part2_3();
    
    std::cout << "\nAll tests completed successfully!\n";
    return 0;
}