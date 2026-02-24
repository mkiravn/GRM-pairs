/*
 * Optimized regression using Cholesky decomposition
 * Much faster than Gaussian elimination for normal equations
 */

/* Fast Cholesky decomposition for positive definite matrices */
static int cholesky_decomp(double* A, uint32_t n) {
    for (uint32_t i = 0; i < n; i++) {
        for (uint32_t j = 0; j < i; j++) {
            double sum = 0.0;
            for (uint32_t k = 0; k < j; k++) {
                sum += A[i * n + k] * A[j * n + k];
            }
            A[i * n + j] = (A[i * n + j] - sum) / A[j * n + j];
        }
        
        double sum = 0.0;
        for (uint32_t k = 0; k < i; k++) {
            sum += A[i * n + k] * A[i * n + k];
        }
        double diag = A[i * n + i] - sum;
        if (diag <= 0.0) return 1; /* Not positive definite */
        A[i * n + i] = sqrt(diag);
    }
    return 0;
}

/* Solve Ly = b where L is lower triangular */
static void forward_solve(const double* L, const double* b, double* y, uint32_t n) {
    for (uint32_t i = 0; i < n; i++) {
        double sum = 0.0;
        for (uint32_t j = 0; j < i; j++) {
            sum += L[i * n + j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i * n + i];
    }
}

/* Solve L^T x = y where L^T is upper triangular */
static void backward_solve(const double* L, const double* y, double* x, uint32_t n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (uint32_t j = i + 1; j < n; j++) {
            sum += L[j * n + i] * x[j]; /* L^T[i][j] = L[j][i] */
        }
        x[i] = (y[i] - sum) / L[i * n + i];
    }
}
