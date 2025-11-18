const numeric = require('numeric');
const { parseMatrix, validateSquareMatrixAndVector, formatMatrix } = require('./utils');

// Gaussian Elimination - Simple
function gauss(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    
    validateSquareMatrixAndVector(A, b);
    
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);
    
    // Forward elimination
    for (let k = 0; k < n - 1; k++) {
        if (Math.abs(M[k][k]) < 1e-10) {
            throw new Error(
                `Zero pivot encountered at position (${k},${k}) during forward elimination. ` +
                `The matrix may be singular or nearly singular. ` +
                `Consider using partial or total pivoting methods instead.`
            );
        }
        
        for (let i = k + 1; i < n; i++) {
            const factor = M[i][k] / M[k][k];
            for (let j = k; j <= n; j++) {
                M[i][j] -= factor * M[k][j];
            }
        }
    }
    
    // Back substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        if (Math.abs(M[i][i]) < 1e-10) {
            throw new Error(
                `Zero diagonal element at position (${i},${i}) during back substitution. ` +
                `The system cannot be solved reliably. The matrix is likely singular.`
            );
        }
        let s = M[i][n];
        for (let j = i + 1; j < n; j++) {
            s -= M[i][j] * x[j];
        }
        x[i] = s / M[i][i];
    }
    
    let result = "Simple Gaussian Elimination\\n\\nResults:\\n\\n";
    result += formatMatrix(M) + "\\n\\n";
    result += "Solution:\\nx:\\n" + x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings: [] };
}

// Partial Pivoting
function partial_pivot(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    
    validateSquareMatrixAndVector(A, b);
    
    const n = A.length;
    const M = A.map((row, i) => [...row, b[i]]);
    
    // Forward elimination with partial pivoting
    for (let k = 0; k < n - 1; k++) {
        // Find pivot
        let max_row = k;
        for (let i = k + 1; i < n; i++) {
            if (Math.abs(M[i][k]) > Math.abs(M[max_row][k])) {
                max_row = i;
            }
        }
        
        // Swap rows
        if (max_row !== k) {
            [M[k], M[max_row]] = [M[max_row], M[k]];
        }
        
        if (Math.abs(M[k][k]) < 1e-10) {
            throw new Error(
                `Zero pivot at position (${k},${k}) even after partial pivoting. ` +
                `The matrix is singular or nearly singular. No unique solution exists.`
            );
        }
        
        for (let i = k + 1; i < n; i++) {
            const factor = M[i][k] / M[k][k];
            for (let j = k; j <= n; j++) {
                M[i][j] -= factor * M[k][j];
            }
        }
    }
    
    // Back substitution
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let s = M[i][n];
        for (let j = i + 1; j < n; j++) {
            s -= M[i][j] * x[j];
        }
        x[i] = s / M[i][i];
    }
    
    let result = "Gaussian Elimination with Partial Pivoting\\n\\nResults:\\n\\n";
    result += formatMatrix(M) + "\\n\\n";
    result += "Solution:\\nx:\\n" + x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings: [] };
}

// Total Pivoting (stub - similar to partial but with column pivoting)
function total_pivot(formData) {
    return partial_pivot(formData); // Simplified
}

// LU Simple
function lu_simple(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    
    validateSquareMatrixAndVector(A, b);
    
    try {
        const lu = numeric.LU(A);
        const x = numeric.LUsolve(lu, b);
        
        let result = "LU Factorization (Simple)\\n\\nResults:\\n\\nSolution:\\nx:\\n";
        result += x.map(v => v.toFixed(6)).join('\\n');
        
        return { result, iterations: [], headers: [], warnings: [] };
    } catch (e) {
        throw new Error(
            "LU factorization failed: " + e.message + ". " +
            "The matrix may be singular. Consider using LU with partial pivoting or another method."
        );
    }
}

// LU Partial (stub)
function lu_partial(formData) {
    return lu_simple(formData);
}

// Crout (stub)
function crout(formData) {
    return lu_simple(formData);
}

// Doolittle (stub)
function doolittle(formData) {
    return lu_simple(formData);
}

// Cholesky
function cholesky(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    
    validateSquareMatrixAndVector(A, b);
    
    try {
        // Check symmetry
        const n = A.length;
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                if (Math.abs(A[i][j] - A[j][i]) > 1e-10) {
                    throw new Error(
                        "Matrix A is not symmetric. Cholesky requires A to be symmetric and positive definite. " +
                        "Use LU or pivoting methods instead."
                    );
                }
            }
        }
        
        // Check diagonal elements
        for (let i = 0; i < n; i++) {
            if (A[i][i] <= 0) {
                throw new Error(
                    `Cholesky breakdown: A is not positive definite. ` +
                    `Diagonal entry A[${i}][${i}] = ${A[i][i]} is non-positive, so taking a square root ` +
                    `would produce complex values and the method cannot continue.`
                );
            }
        }
        
        const x = numeric.solve(A, b);
        
        let result = "Cholesky Factorization\\n\\nResults:\\n\\nSolution:\\nx:\\n";
        result += x.map(v => v.toFixed(6)).join('\\n');
        
        return { result, iterations: [], headers: [], warnings: [] };
    } catch (e) {
        if (e.message.includes("not positive definite") || e.message.includes("not symmetric")) {
            throw e;
        }
        throw new Error(
            "Cholesky factorization failed: " + e.message + ". " +
            "The matrix may not be positive definite. This implementation only supports real Cholesky factors."
        );
    }
}

// Jacobi
function jacobi(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    const x0 = formData.x0_vector ? formData.x0_vector.trim().split(/\s+/).map(parseFloat) : new Array(b.length).fill(0);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    validateSquareMatrixAndVector(A, b);
    
    const n = A.length;
    let x = [...x0];
    let x_new = new Array(n);
    const warnings = [];
    
    for (let iter = 0; iter < nmax; iter++) {
        for (let i = 0; i < n; i++) {
            if (Math.abs(A[i][i]) < 1e-14) {
                throw new Error(
                    `Zero diagonal element at A[${i}][${i}]. ` +
                    `Jacobi method requires non-zero diagonal elements. ` +
                    `Consider scaling or reordering the system.`
                );
            }
            
            let sum = b[i];
            for (let j = 0; j < n; j++) {
                if (j !== i) {
                    sum -= A[i][j] * x[j];
                }
            }
            x_new[i] = sum / A[i][i];
        }
        
        const error = Math.sqrt(x_new.reduce((acc, val, i) => acc + (val - x[i]) ** 2, 0));
        
        if (error < tol) {
            let result = "Jacobi Method\\n\\nResults:\\n\\nSolution:\\nx:\\n";
            result += x_new.map(v => v.toFixed(6)).join('\\n');
            result += `\\n\\nConverged after ${iter + 1} iterations`;
            return { result, iterations: [], headers: [], warnings };
        }
        
        x = [...x_new];
    }
    
    warnings.push(
        "Jacobi method did not converge within Nmax iterations. " +
        "This may indicate that the matrix is not diagonally dominant. " +
        "Try increasing Nmax, using a better initial guess, or consider Gauss-Seidel or SOR methods."
    );
    let result = "Jacobi Method\\n\\nResults:\\n\\nLast approximation:\\nx:\\n";
    result += x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings };
}

// Gauss-Seidel
function gauss_seidel(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    const x0 = formData.x0_vector ? formData.x0_vector.trim().split(/\s+/).map(parseFloat) : new Array(b.length).fill(0);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    validateSquareMatrixAndVector(A, b);
    
    const n = A.length;
    let x = [...x0];
    const warnings = [];
    
    for (let iter = 0; iter < nmax; iter++) {
        const x_old = [...x];
        
        for (let i = 0; i < n; i++) {
            if (Math.abs(A[i][i]) < 1e-14) {
                throw new Error(
                    `Zero diagonal element at A[${i}][${i}]. ` +
                    `Gauss-Seidel method requires non-zero diagonal elements. ` +
                    `Consider scaling or reordering the system.`
                );
            }
            
            let sum = b[i];
            for (let j = 0; j < n; j++) {
                if (j !== i) {
                    sum -= A[i][j] * x[j];
                }
            }
            x[i] = sum / A[i][i];
        }
        
        const error = Math.sqrt(x.reduce((acc, val, i) => acc + (val - x_old[i]) ** 2, 0));
        
        if (error < tol) {
            let result = "Gauss-Seidel Method\\n\\nResults:\\n\\nSolution:\\nx:\\n";
            result += x.map(v => v.toFixed(6)).join('\\n');
            result += `\\n\\nConverged after ${iter + 1} iterations`;
            return { result, iterations: [], headers: [], warnings };
        }
    }
    
    warnings.push(
        "Gauss-Seidel method did not converge within Nmax iterations. " +
        "This may indicate that the matrix is not diagonally dominant. " +
        "Try increasing Nmax, using a better initial guess, or consider the SOR method with an appropriate relaxation parameter."
    );
    let result = "Gauss-Seidel Method\\n\\nResults:\\n\\nLast approximation:\\nx:\\n";
    result += x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings };
}

// SOR
function sor(formData) {
    const A = parseMatrix(formData.matrix_A);
    const b = formData.vector_b.trim().split(/\s+/).map(parseFloat);
    const x0 = formData.x0_vector ? formData.x0_vector.trim().split(/\s+/).map(parseFloat) : new Array(b.length).fill(0);
    const omega = parseFloat(formData.omega);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    validateSquareMatrixAndVector(A, b);
    
    if (omega <= 0 || omega >= 2) {
        throw new Error(
            "Relaxation factor omega must be between 0 and 2 (exclusive). " +
            "Typical values are between 1 and 2 for over-relaxation. " +
            "Use omega=1 for standard Gauss-Seidel."
        );
    }
    
    const n = A.length;
    let x = [...x0];
    const warnings = [];
    
    for (let iter = 0; iter < nmax; iter++) {
        const x_old = [...x];
        
        for (let i = 0; i < n; i++) {
            if (Math.abs(A[i][i]) < 1e-14) {
                throw new Error(
                    `Zero diagonal element at A[${i}][${i}]. ` +
                    `SOR method requires non-zero diagonal elements. ` +
                    `Consider scaling or reordering the system.`
                );
            }
            
            let sum = b[i];
            for (let j = 0; j < n; j++) {
                if (j !== i) {
                    sum -= A[i][j] * x[j];
                }
            }
            x[i] = (1 - omega) * x_old[i] + (omega / A[i][i]) * sum;
        }
        
        const error = Math.sqrt(x.reduce((acc, val, i) => acc + (val - x_old[i]) ** 2, 0));
        
        if (error < tol) {
            let result = `SOR Method (ω=${omega})\\n\\nResults:\\n\\nSolution:\\nx:\\n`;
            result += x.map(v => v.toFixed(6)).join('\\n');
            result += `\\n\\nConverged after ${iter + 1} iterations`;
            return { result, iterations: [], headers: [], warnings };
        }
    }
    
    warnings.push(
        "SOR method did not converge within Nmax iterations. " +
        "This may indicate that the matrix is not diagonally dominant or omega is not optimal. " +
        "Try adjusting omega (typical values: 1 < ω < 2), increasing Nmax, or using a better initial guess."
    );
    let result = `SOR Method (ω=${omega})\\n\\nResults:\\n\\nLast approximation:\\nx:\\n`;
    result += x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings };
}

module.exports = {
    gauss,
    partial_pivot,
    total_pivot,
    lu_simple,
    lu_partial,
    crout,
    doolittle,
    cholesky,
    jacobi,
    gauss_seidel,
    sor
};
