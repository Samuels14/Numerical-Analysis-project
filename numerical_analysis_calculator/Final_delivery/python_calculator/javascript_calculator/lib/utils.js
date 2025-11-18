const math = require('mathjs');

// Create a safe function evaluator from expression string
function createFunction(expr) {
    // Normalize expression for mathjs:
    // 1. Replace ln with log (mathjs uses log for natural logarithm)
    expr = expr.replace(/ln\(/g, 'log(');
    expr = expr.replace(/Log\(/g, 'log(');
    
    // 2. Replace ** with ^ for exponentiation (Python syntax to mathjs)
    expr = expr.replace(/\*\*/g, '^');
    
    // 3. Replace common Python log functions
    // log10(x) works in mathjs, but ensure no conflicts
    expr = expr.replace(/log10\(/g, 'log10(');
    
    try {
        const compiled = math.compile(expr);
        return (x) => {
            try {
                return compiled.evaluate({ x });
            } catch (e) {
                throw new Error(`Error evaluating function at x = ${x}: ${e.message}`);
            }
        };
    } catch (e) {
        throw new Error(`Invalid expression: ${e.message}`);
    }
}

// Safe function evaluation with validation
function safeEval(f, x, context = "f") {
    try {
        const y = f(x);
        if (isNaN(y) || !isFinite(y)) {
            throw new Error(`Evaluation of ${context}(x) at x = ${x} produced an invalid value (NaN/Inf).`);
        }
        return y;
    } catch (e) {
        throw new Error(`Error evaluating ${context}(x) at x = ${x}: ${e.message}`);
    }
}

// Numeric derivative using central differences
function numericDerivative(f, x, h = 1e-6) {
    return (safeEval(f, x + h) - safeEval(f, x - h)) / (2 * h);
}

// Numeric second derivative
function numericSecondDerivative(f, x, h = 1e-5) {
    return (safeEval(f, x + h) - 2 * safeEval(f, x) + safeEval(f, x - h)) / (h * h);
}

// Format a matrix for display
function formatMatrix(matrix) {
    return matrix.map(row => 
        ' ' + row.map(x => x.toFixed(6)).join('  ') + ' '
    ).join('\n');
}

// Parse matrix from text (rows separated by newlines, values by spaces)
function parseMatrix(matrixText) {
    try {
        const rows = matrixText.trim().split('\n').filter(row => row.trim());
        const matrix = rows.map(row => 
            row.trim().split(/\s+/).map(x => parseFloat(x))
        );
        
        if (matrix.length === 0) {
            throw new Error("Matrix cannot be empty.");
        }
        
        const numCols = matrix[0].length;
        for (const row of matrix) {
            if (row.length !== numCols) {
                throw new Error("All rows of the matrix must have the same number of columns.");
            }
        }
        
        return matrix;
    } catch (e) {
        throw new Error("Error reading matrix. Use only numbers separated by spaces, one row per line.");
    }
}

// Validate square matrix and vector
function validateSquareMatrixAndVector(A, b) {
    const n = A.length;
    if (n === 0) {
        throw new Error("Matrix A cannot be empty.");
    }
    for (const row of A) {
        if (row.length !== n) {
            throw new Error("Matrix A must be square (same number of rows and columns).");
        }
    }
    if (b.length !== n) {
        throw new Error("Vector b must have the same length as the number of rows in A.");
    }
}

module.exports = {
    createFunction,
    safeEval,
    numericDerivative,
    numericSecondDerivative,
    formatMatrix,
    parseMatrix,
    validateSquareMatrixAndVector
};
