const { parseMatrix } = require('./utils');

// Forward Substitution
function forward_substitution(formData) {
    const M = parseMatrix(formData.augmented_matrix);
    const n = M.length;
    const x = new Array(n).fill(0);
    
    if (Math.abs(M[0][0]) < 1e-14) {
        throw new Error("Zero diagonal element in L[0][0].");
    }
    
    x[0] = M[0][n] / M[0][0];
    
    for (let i = 1; i < n; i++) {
        if (Math.abs(M[i][i]) < 1e-14) {
            throw new Error(`Zero diagonal element in L[${i}][${i}].`);
        }
        let s = M[i][n];
        for (let j = 0; j < i; j++) {
            s -= M[i][j] * x[j];
        }
        x[i] = s / M[i][i];
    }
    
    let result = "Forward Substitution\\n\\nResults:\\n\\nSolution:\\nx:\\n";
    result += x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings: [] };
}

// Backward Substitution
function backward_substitution(formData) {
    const M = parseMatrix(formData.augmented_matrix);
    const n = M.length;
    const x = new Array(n).fill(0);
    
    if (Math.abs(M[n-1][n-1]) < 1e-14) {
        throw new Error(`Zero diagonal element in U[${n-1}][${n-1}].`);
    }
    
    x[n-1] = M[n-1][n] / M[n-1][n-1];
    
    for (let i = n-2; i >= 0; i--) {
        if (Math.abs(M[i][i]) < 1e-14) {
            throw new Error(`Zero diagonal element in U[${i}][${i}].`);
        }
        let s = M[i][n];
        for (let j = i+1; j < n; j++) {
            s -= M[i][j] * x[j];
        }
        x[i] = s / M[i][i];
    }
    
    let result = "Backward Substitution\\n\\nResults:\\n\\nSolution:\\nx:\\n";
    result += x.map(v => v.toFixed(6)).join('\\n');
    
    return { result, iterations: [], headers: [], warnings: [] };
}

module.exports = {
    forward_substitution,
    backward_substitution
};
