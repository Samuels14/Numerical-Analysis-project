const numeric = require('numeric');

// Helper function to format matrix for display
function formatMatrix(matrix) {
    return matrix.map(row => 
        ' ' + row.map(x => x.toFixed(6)).join('  ') + ' '
    ).join('\n');
}

// Vandermonde
function vandermonde(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 2) {
        throw new Error("At least 2 points are required for interpolation.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const A = [];
    
    // Build Vandermonde matrix
    for (let i = 0; i < n; i++) {
        const row = [];
        for (let j = 0; j < n; j++) {
            row.push(Math.pow(X[i], n - 1 - j));
        }
        A.push(row);
    }
    
    if (n > 10) {
        warnings.push("Using many points in a Vandermonde system can produce numerical instability (Runge phenomenon).");
    }
    
    const coef = numeric.solve(A, Y);
    const degree = n - 1;
    
    let result = "Vandermonde\\n\\nResults:\\n\\nVandermonde matrix:\\n\\n";
    result += formatMatrix(A) + "\\n";
    result += "\\nPolynomial coefficients:\\n\\n";
    
    for (const c of coef) {
        result += ` ${c.toFixed(6)} `;
    }
    
    result += "\\n\\nPolynomial:\\n\\n";
    
    const terms = coef.map((c, i) => {
        const exp = degree - i;
        if (exp > 1) return `${c.toFixed(6)}x^${exp}`;
        if (exp === 1) return `${c.toFixed(6)}x`;
        return `${c.toFixed(6)}`;
    });
    
    result += terms.join('+').replace(/\+-/g, '-');
    
    return { result, iterations: [], headers: [], warnings };
}

// Divided Differences (Newton)
function divided_differences(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 2) {
        throw new Error("At least 2 points are required for interpolation.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const D = Array(n).fill(0).map(() => Array(n).fill(0));
    
    // First column
    for (let i = 0; i < n; i++) {
        D[i][0] = Y[i];
    }
    
    // Compute divided differences
    for (let j = 1; j < n; j++) {
        for (let i = j; i < n; i++) {
            D[i][j] = (D[i][j-1] - D[i-1][j-1]) / (X[i] - X[i-j]);
        }
    }
    
    const coef = D.map((row, i) => row[i]);
    const degree = n - 1;
    
    let result = "Newton\\n\\nResults:\\n\\nDivided differences table:\\n\\n";
    result += formatMatrix(D) + "\\n";
    result += "\\nNewton polynomial coefficients:\\n\\n";
    
    for (const c of coef) {
        result += ` ${c.toFixed(6)} `;
    }
    
    result += "\\n\\nNewton polynomial:\\n\\n";
    
    const parts = [`${coef[0].toFixed(6)}`];
    for (let i = 1; i < coef.length; i++) {
        const factors = X.slice(0, i).map(x => 
            x !== 0 ? `(x-(${x.toFixed(6)}))` : 'x'
        ).join('');
        parts.push(`${coef[i].toFixed(6)}${factors}`);
    }
    
    result += parts.join('+').replace(/\+-/g, '-').replace(/\(x-\(-/g, '(x+').replace(/\)\)/g, ')');
    
    return { result, iterations: [], headers: [], warnings };
}

// Lagrange
function lagrange(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 2) {
        throw new Error("At least 2 points are required for interpolation.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const L = Array(n).fill(0).map(() => Array(n).fill(0));
    
    // Compute Lagrange basis polynomials
    for (let i = 0; i < n; i++) {
        const aux0 = X.filter((_, j) => j !== i);
        let aux = [1.0, -aux0[0]];
        
        for (let j = 1; j < n - 1; j++) {
            // Polynomial multiplication (convolution)
            const newAux = Array(aux.length + 1).fill(0);
            for (let k = 0; k < aux.length; k++) {
                newAux[k] += aux[k];
                newAux[k + 1] -= aux[k] * aux0[j];
            }
            aux = newAux;
        }
        
        // Evaluate polynomial at X[i] for denominator
        let denominator = 0;
        for (let k = 0; k < aux.length; k++) {
            denominator += aux[k] * Math.pow(X[i], aux.length - 1 - k);
        }
        
        // Normalize
        for (let k = 0; k < aux.length; k++) {
            L[i][k] = aux[k] / denominator;
        }
        
        // Pad with zeros if needed
        while (L[i].length < n) {
            L[i].push(0);
        }
    }
    
    const degree = n - 1;
    
    let result = "Lagrange\\n\\nResults:\\n\\nLagrange interpolating polynomials:\\n\\n";
    
    for (let i = 0; i < L.length; i++) {
        const Li = L[i];
        const Li_degree = Li.length - 1;
        const terms = [];
        
        for (let j = 0; j < Li.length; j++) {
            const exp = Li_degree - j;
            if (Math.abs(Li[j]) > 1e-10) {  // Only show non-zero terms
                if (exp > 1) {
                    terms.push(`${Li[j].toFixed(6)}x^${exp}`);
                } else if (exp === 1) {
                    terms.push(`${Li[j].toFixed(6)}x`);
                } else {
                    terms.push(`${Li[j].toFixed(6)}`);
                }
            }
        }
        
        result += terms.join('+').replace(/\+-/g, '-') + `   //L${i}\\n`;
    }
    
    result += "\\n\\nPolynomial:\\n\\n";
    const poly_parts = Y.map((y, i) => `${y.toFixed(1)}*L${i}`);
    result += poly_parts.join('+').replace(/\+-/g, '-');
    
    return { result, iterations: [], headers: [], warnings };
}

// Linear Spline
function linear_spline(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 2) {
        throw new Error("At least 2 points are required.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const A = [];
    const B = [];
    
    // Calculate coefficients for each segment
    for (let i = 0; i < n - 1; i++) {
        const a = Y[i];
        const b = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
        A.push(a);
        B.push(b);
    }
    
    let result = "Linear Spline\\n\\nResults:\\n\\nCoefficients:\\n\\n";
    
    // Create coefficient table
    result += "Interval | a_i | b_i\\n";
    result += "---------|-----|-----\\n";
    
    for (let i = 0; i < A.length; i++) {
        result += `${i} | ${A[i].toFixed(6)} | ${B[i].toFixed(6)}\\n`;
    }
    
    result += "\\n\\nSpline equations:\\n\\n";
    
    for (let i = 0; i < A.length; i++) {
        const a_str = A[i].toFixed(6);
        const b_str = B[i] >= 0 ? `+${B[i].toFixed(6)}` : `${B[i].toFixed(6)}`;
        const x_val = X[i].toFixed(6);
        result += `S${i}(x) = ${a_str} ${b_str}(x-${x_val})   for x in [${X[i].toFixed(2)}, ${X[i+1].toFixed(2)}]\\n`;
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Quadratic Spline (simplified stub)
function quadratic_spline(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 3) {
        throw new Error("At least 3 points are required for quadratic splines.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const m = 3 * (n - 1);
    
    // Build system of equations for quadratic splines
    const A = Array(m).fill(0).map(() => Array(m).fill(0));
    const b = Array(m).fill(0);
    
    // Equations: S_i(x_{i+1}) = y_{i+1}
    for (let i = 0; i < n - 1; i++) {
        const row_idx = i + 1;
        const col_start = 3 * i;
        A[row_idx][col_start] = X[i + 1] ** 2;
        A[row_idx][col_start + 1] = X[i + 1];
        A[row_idx][col_start + 2] = 1;
        b[row_idx] = Y[i + 1];
    }
    
    // S_0(x_0) = y_0
    A[0][0] = X[0] ** 2;
    A[0][1] = X[0];
    A[0][2] = 1;
    b[0] = Y[0];
    
    // Continuity at interior points: S_{i-1}(x_i) = S_i(x_i)
    for (let i = 1; i < n - 1; i++) {
        const row_idx = n - 1 + i;
        const col_start = 3 * i - 3;
        A[row_idx][col_start] = X[i] ** 2;
        A[row_idx][col_start + 1] = X[i];
        A[row_idx][col_start + 2] = 1;
        A[row_idx][col_start + 3] = -(X[i] ** 2);
        A[row_idx][col_start + 4] = -X[i];
        A[row_idx][col_start + 5] = -1;
        b[row_idx] = 0;
    }
    
    // First derivative continuity: S'_{i-1}(x_i) = S'_i(x_i)
    for (let i = 1; i < n - 1; i++) {
        const row_idx = 2 * n - 3 + i;
        const col_start = 3 * i - 3;
        A[row_idx][col_start] = 2 * X[i];
        A[row_idx][col_start + 1] = 1;
        A[row_idx][col_start + 3] = -2 * X[i];
        A[row_idx][col_start + 4] = -1;
        b[row_idx] = 0;
    }
    
    // S''_0(x_0) = 0 (natural spline boundary condition)
    A[m - 1][0] = 2;
    b[m - 1] = 0;
    
    // Solve using Gaussian elimination with partial pivoting
    const numeric = require('numeric');
    let Saux;
    try {
        Saux = numeric.solve(A, b);
    } catch (e) {
        throw new Error("Failed to solve the system for quadratic spline coefficients.");
    }
    
    const Coef = [];
    for (let i = 0; i < n - 1; i++) {
        Coef.push([Saux[3 * i], Saux[3 * i + 1], Saux[3 * i + 2]]);
    }
    
    let result = "Quadratic Splines\\n\\nResults:\\n\\nSpline coefficients:\\n\\n";
    
    // Create coefficient table
    result += "Interval | a_i | b_i | c_i\\n";
    result += "---------|-----|-----|-----\\n";
    
    for (let i = 0; i < Coef.length; i++) {
        const [a_i, b_i, c_i] = Coef[i];
        result += `${i} | ${a_i.toFixed(6)} | ${b_i.toFixed(6)} | ${c_i.toFixed(6)}\\n`;
    }
    
    result += "\\n\\nSplines:\\n\\n";
    
    for (let i = 0; i < Coef.length; i++) {
        const [a_i, b_i, c_i] = Coef[i];
        const b_sign = b_i >= 0 ? "+" : "";
        const c_sign = c_i >= 0 ? "+" : "";
        result += `S${i}(x) = ${a_i.toFixed(6)}x^2${b_sign}${b_i.toFixed(6)}x${c_sign}${c_i.toFixed(6)}\\n`;
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Cubic Spline (simplified stub)
function cubic_spline(formData) {
    const X = formData.x_values.trim().split(/\s+/).map(parseFloat);
    const Y = formData.y_values.trim().split(/\s+/).map(parseFloat);
    
    const warnings = [];
    
    if (X.length !== Y.length) {
        throw new Error("X and Y must have the same number of points.");
    }
    if (X.length < 3) {
        throw new Error("At least 3 points are required for cubic splines.");
    }
    
    // Check for duplicate X values
    const xSet = new Set(X);
    if (xSet.size < X.length) {
        throw new Error("Duplicate X values found. All X values must be unique for interpolation methods.");
    }
    
    // Check if X values are in ascending order
    for (let i = 0; i < X.length - 1; i++) {
        if (X[i] >= X[i + 1]) {
            throw new Error("X values must be entered in strictly ascending order (from smallest to largest).");
        }
    }
    
    const n = X.length;
    const m = 4 * (n - 1);
    
    // Build system of equations for cubic splines
    const A = Array(m).fill(0).map(() => Array(m).fill(0));
    const b = Array(m).fill(0);
    
    // Equations: S_i(x_{i+1}) = y_{i+1}
    for (let i = 0; i < n - 1; i++) {
        const row_idx = i + 1;
        const col_start = 4 * i;
        A[row_idx][col_start] = X[i + 1] ** 3;
        A[row_idx][col_start + 1] = X[i + 1] ** 2;
        A[row_idx][col_start + 2] = X[i + 1];
        A[row_idx][col_start + 3] = 1;
        b[row_idx] = Y[i + 1];
    }
    
    // S_0(x_0) = y_0
    A[0][0] = X[0] ** 3;
    A[0][1] = X[0] ** 2;
    A[0][2] = X[0];
    A[0][3] = 1;
    b[0] = Y[0];
    
    // Continuity at interior points: S_{i-1}(x_i) = S_i(x_i)
    for (let i = 1; i < n - 1; i++) {
        const row_idx = n - 1 + i;
        const col_start = 4 * i - 4;
        A[row_idx][col_start] = X[i] ** 3;
        A[row_idx][col_start + 1] = X[i] ** 2;
        A[row_idx][col_start + 2] = X[i];
        A[row_idx][col_start + 3] = 1;
        A[row_idx][col_start + 4] = -(X[i] ** 3);
        A[row_idx][col_start + 5] = -(X[i] ** 2);
        A[row_idx][col_start + 6] = -X[i];
        A[row_idx][col_start + 7] = -1;
        b[row_idx] = 0;
    }
    
    // First derivative continuity: S'_{i-1}(x_i) = S'_i(x_i)
    for (let i = 1; i < n - 1; i++) {
        const row_idx = 2 * n - 3 + i;
        const col_start = 4 * i - 4;
        A[row_idx][col_start] = 3 * X[i] ** 2;
        A[row_idx][col_start + 1] = 2 * X[i];
        A[row_idx][col_start + 2] = 1;
        A[row_idx][col_start + 4] = -3 * X[i] ** 2;
        A[row_idx][col_start + 5] = -2 * X[i];
        A[row_idx][col_start + 6] = -1;
        b[row_idx] = 0;
    }
    
    // Second derivative continuity: S''_{i-1}(x_i) = S''_i(x_i)
    for (let i = 1; i < n - 1; i++) {
        const row_idx = 3 * n - 5 + i;
        const col_start = 4 * i - 4;
        A[row_idx][col_start] = 6 * X[i];
        A[row_idx][col_start + 1] = 2;
        A[row_idx][col_start + 4] = -6 * X[i];
        A[row_idx][col_start + 5] = -2;
        b[row_idx] = 0;
    }
    
    // Natural spline boundary conditions: S''_0(x_0) = 0
    A[m - 2][0] = 6 * X[0];
    A[m - 2][1] = 2;
    b[m - 2] = 0;
    
    // S''_{n-2}(x_{n-1}) = 0
    A[m - 1][m - 4] = 6 * X[n - 1];
    A[m - 1][m - 3] = 2;
    b[m - 1] = 0;
    
    // Solve using Gaussian elimination with partial pivoting
    const numeric = require('numeric');
    let Saux;
    try {
        Saux = numeric.solve(A, b);
    } catch (e) {
        throw new Error("Failed to solve the system for cubic spline coefficients.");
    }
    
    const Coef = [];
    for (let i = 0; i < n - 1; i++) {
        Coef.push([Saux[4 * i], Saux[4 * i + 1], Saux[4 * i + 2], Saux[4 * i + 3]]);
    }
    
    let result = "Cubic Splines\\n\\nResults:\\n\\nSpline coefficients:\\n\\n";
    
    // Create coefficient table
    result += "Interval | a_i | b_i | c_i | d_i\\n";
    result += "---------|-----|-----|-----|-----\\n";
    
    for (let i = 0; i < Coef.length; i++) {
        const [a_i, b_i, c_i, d_i] = Coef[i];
        result += `${i} | ${a_i.toFixed(6)} | ${b_i.toFixed(6)} | ${c_i.toFixed(6)} | ${d_i.toFixed(6)}\\n`;
    }
    
    result += "\\n\\nSplines:\\n\\n";
    
    for (let i = 0; i < Coef.length; i++) {
        const [a_i, b_i, c_i, d_i] = Coef[i];
        const b_sign = b_i >= 0 ? "+" : "";
        const c_sign = c_i >= 0 ? "+" : "";
        const d_sign = d_i >= 0 ? "+" : "";
        result += `S${i}(x) = ${a_i.toFixed(6)}x^3${b_sign}${b_i.toFixed(6)}x^2${c_sign}${c_i.toFixed(6)}x${d_sign}${d_i.toFixed(6)}\\n`;
    }
    
    return { result, iterations: [], headers: [], warnings };
}

module.exports = {
    vandermonde,
    divided_differences,
    lagrange,
    linear_spline,
    quadratic_spline,
    cubic_spline
};
