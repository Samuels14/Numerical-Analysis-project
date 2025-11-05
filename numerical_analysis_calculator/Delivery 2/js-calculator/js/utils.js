/**
 * UTILIDADES PARA CÁLCULOS NUMÉRICOS
 * Equivalente a las funciones auxiliares de Python
 */

// Configuración global
const CONFIG = {
  EPSILON: 1e-10,
  MAX_ITERATIONS: 1000,
  DEFAULT_TOL: 1e-7
};

/**
 * Crea una función evaluable desde una expresión de texto
 * Equivalente a create_function en Python
 * @param {string} expr - Expresión matemática
 * @returns {Function} - Función evaluable
 */
function createFunction(expr) {
  try {
    // Reemplazar ln y Log por log (natural)
    expr = expr.replace(/ln\(/gi, 'log(');
    expr = expr.replace(/Log\(/gi, 'log(');
    
    // Parsear con math.js
    const node = math.parse(expr);
    
    // Compilar a función
    const compiled = node.compile();
    
    // Retornar función que evalúa en x
    return (x) => {
      try {
        return compiled.evaluate({ x: x });
      } catch (e) {
        throw new Error(`Error evaluating function at x=${x}: ${e.message}`);
      }
    };
  } catch (e) {
    throw new Error(`Invalid expression: ${expr}. Error: ${e.message}`);
  }
}

/**
 * Derivada numérica usando diferencias centradas
 * @param {Function} f - Función a derivar
 * @param {number} x - Punto de evaluación
 * @param {number} h - Tamaño de paso
 * @returns {number} - Valor de la derivada
 */
function numericDerivative(f, x, h = 1e-6) {
  return (f(x + h) - f(x - h)) / (2 * h);
}

/**
 * Segunda derivada numérica
 * @param {Function} f - Función a derivar
 * @param {number} x - Punto de evaluación
 * @param {number} h - Tamaño de paso
 * @returns {number} - Valor de la segunda derivada
 */
function numericSecondDerivative(f, x, h = 1e-5) {
  return (f(x + h) - 2 * f(x) + f(x - h)) / (h ** 2);
}

/**
 * Parsea una matriz desde texto
 * @param {string} matrixText - Texto con la matriz
 * @returns {Array<Array<number>>} - Matriz como array 2D
 */
function parseMatrix(matrixText) {
  const rows = matrixText.trim().split('\n').filter(row => row.trim());
  return rows.map(row => 
    row.trim().split(/\s+/).map(x => parseFloat(x))
  );
}

/**
 * Parsea un vector desde texto
 * @param {string} vectorText - Texto con el vector
 * @returns {Array<number>} - Vector como array
 */
function parseVector(vectorText) {
  return vectorText.trim().split(/\s+/).map(x => parseFloat(x));
}

/**
 * Formatea una matriz para mostrar
 * @param {Array<Array<number>>} matrix - Matriz a formatear
 * @param {number} decimals - Número de decimales
 * @returns {string} - Matriz formateada
 */
function formatMatrix(matrix, decimals = 6) {
  return matrix.map(row => 
    ' ' + row.map(x => x.toFixed(decimals)).join('  ') + ' '
  ).join('\n');
}

/**
 * Formatea un número en notación científica si es muy pequeño
 * @param {number} num - Número a formatear
 * @param {number} decimals - Decimales
 * @returns {string} - Número formateado
 */
function formatNumber(num, decimals = 6) {
  if (Math.abs(num) < 1e-3 && num !== 0) {
    return num.toExponential(1);
  }
  return num.toFixed(decimals);
}

/**
 * Copia profunda de una matriz
 * @param {Array<Array<number>>} matrix - Matriz a copiar
 * @returns {Array<Array<number>>} - Copia de la matriz
 */
function copyMatrix(matrix) {
  return matrix.map(row => [...row]);
}

/**
 * Crea una matriz identidad
 * @param {number} n - Tamaño de la matriz
 * @returns {Array<Array<number>>} - Matriz identidad
 */
function identityMatrix(n) {
  const I = [];
  for (let i = 0; i < n; i++) {
    I[i] = [];
    for (let j = 0; j < n; j++) {
      I[i][j] = (i === j) ? 1.0 : 0.0;
    }
  }
  return I;
}

/**
 * Crea una matriz de ceros
 * @param {number} rows - Número de filas
 * @param {number} cols - Número de columnas
 * @returns {Array<Array<number>>} - Matriz de ceros
 */
function zerosMatrix(rows, cols) {
  const M = [];
  for (let i = 0; i < rows; i++) {
    M[i] = [];
    for (let j = 0; j < cols; j++) {
      M[i][j] = 0.0;
    }
  }
  return M;
}

/**
 * Multiplica matriz por vector
 * @param {Array<Array<number>>} A - Matriz
 * @param {Array<number>} x - Vector
 * @returns {Array<number>} - Resultado
 */
function matrixVectorMultiply(A, x) {
  const n = A.length;
  const result = new Array(n).fill(0);
  
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result[i] += A[i][j] * x[j];
    }
  }
  
  return result;
}

/**
 * Norma euclidiana de un vector
 * @param {Array<number>} v - Vector
 * @returns {number} - Norma
 */
function vectorNorm(v) {
  return Math.sqrt(v.reduce((sum, val) => sum + val * val, 0));
}

/**
 * Resta de dos vectores
 * @param {Array<number>} v1 - Vector 1
 * @param {Array<number>} v2 - Vector 2
 * @returns {Array<number>} - v1 - v2
 */
function vectorSubtract(v1, v2) {
  return v1.map((val, i) => val - v2[i]);
}

/**
 * Valida que una expresión sea válida
 * @param {string} expr - Expresión a validar
 * @returns {boolean} - true si es válida
 */
function validateExpression(expr) {
  try {
    math.parse(expr);
    return true;
  } catch (e) {
    return false;
  }
}

/**
 * Obtiene la diagonal de una matriz
 * @param {Array<Array<number>>} A - Matriz
 * @returns {Array<Array<number>>} - Matriz diagonal
 */
function getDiagonal(A) {
  const n = A.length;
  const D = zerosMatrix(n, n);
  for (let i = 0; i < n; i++) {
    D[i][i] = A[i][i];
  }
  return D;
}

/**
 * Obtiene la parte triangular inferior de una matriz (sin diagonal)
 * @param {Array<Array<number>>} A - Matriz
 * @returns {Array<Array<number>>} - Matriz triangular inferior negada
 */
function getLowerTriangular(A) {
  const n = A.length;
  const L = zerosMatrix(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < i; j++) {
      L[i][j] = -A[i][j];
    }
  }
  return L;
}

/**
 * Obtiene la parte triangular superior de una matriz (sin diagonal)
 * @param {Array<Array<number>>} A - Matriz
 * @returns {Array<Array<number>>} - Matriz triangular superior negada
 */
function getUpperTriangular(A) {
  const n = A.length;
  const U = zerosMatrix(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      U[i][j] = -A[i][j];
    }
  }
  return U;
}

/**
 * Suma de dos matrices
 * @param {Array<Array<number>>} A - Matriz 1
 * @param {Array<Array<number>>} B - Matriz 2
 * @returns {Array<Array<number>>} - A + B
 */
function matrixAdd(A, B) {
  const n = A.length;
  const C = zerosMatrix(n, n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      C[i][j] = A[i][j] + B[i][j];
    }
  }
  return C;
}

/**
 * Inversión de matriz usando math.js
 * @param {Array<Array<number>>} A - Matriz a invertir
 * @returns {Array<Array<number>>} - Matriz inversa
 */
function matrixInverse(A) {
  try {
    return math.inv(A);
  } catch (e) {
    throw new Error('Matrix is singular and cannot be inverted');
  }
}

/**
 * Multiplicación de matrices usando math.js
 * @param {Array<Array<number>>} A - Matriz 1
 * @param {Array<Array<number>>} B - Matriz 2
 * @returns {Array<Array<number>>} - A * B
 */
function matrixMultiply(A, B) {
  return math.multiply(A, B);
}

/**
 * Calcula el radio espectral usando el método de las potencias
 * (método iterativo que no requiere math.eigs)
 * @param {Array<Array<number>>} A - Matriz cuadrada
 * @param {number} maxIter - Iteraciones máximas
 * @returns {number} - Aproximación del radio espectral
 */
function spectralRadiusPowerMethod(A, maxIter = 100) {
  const n = A.length;
  
  // Vector inicial aleatorio
  let v = Array(n).fill(0).map(() => Math.random());
  
  // Normalizar
  let norm = Math.sqrt(v.reduce((sum, val) => sum + val * val, 0));
  v = v.map(val => val / norm);
  
  let lambda = 0;
  
  for (let i = 0; i < maxIter; i++) {
    // v_new = A * v
    const v_new = A.map(row => 
      row.reduce((sum, val, j) => sum + val * v[j], 0)
    );
    
    // Calcular eigenvalor (producto punto v_new · v)
    lambda = v_new.reduce((sum, val, j) => sum + val * v[j], 0);
    
    // Normalizar v_new
    norm = Math.sqrt(v_new.reduce((sum, val) => sum + val * val, 0));
    
    if (norm < 1e-10) break;
    
    v = v_new.map(val => val / norm);
  }
  
  return Math.abs(lambda);
}

/**
 * Calcula el radio espectral (máximo eigenvalor en valor absoluto)
 * Usa método alternativo si math.eigs falla
 * @param {Array<Array<number>>} A - Matriz
 * @returns {number} - Radio espectral
 */
function spectralRadius(A) {
  try {
    // Intentar con math.eigs primero
    const result = math.eigs(A);
    
    // Extraer valores propios
    let eigenvalues;
    
    if (result && result.values) {
      eigenvalues = result.values;
      
      // Si es un objeto Matrix, convertir a array
      if (eigenvalues._data) {
        eigenvalues = eigenvalues._data;
      } else if (typeof eigenvalues.toArray === 'function') {
        eigenvalues = eigenvalues.toArray();
      }
      
      // Si es matriz columna, aplanar
      if (Array.isArray(eigenvalues) && eigenvalues.length > 0 && Array.isArray(eigenvalues[0])) {
        eigenvalues = eigenvalues.map(row => row[0]);
      }
      
      // Calcular magnitudes
      const magnitudes = eigenvalues.map(val => {
        if (typeof val === 'object' && val.re !== undefined) {
          return Math.sqrt(val.re * val.re + (val.im || 0) * (val.im || 0));
        }
        return Math.abs(val);
      });
      
      return Math.max(...magnitudes);
    }
  } catch (e) {
    console.warn('math.eigs failed, using power method:', e.message);
  }
  
  // Si math.eigs falla, usar método de las potencias
  return spectralRadiusPowerMethod(A);
}

// Exportar funciones
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    CONFIG,
    createFunction,
    numericDerivative,
    numericSecondDerivative,
    parseMatrix,
    parseVector,
    formatMatrix,
    formatNumber,
    copyMatrix,
    identityMatrix,
    zerosMatrix,
    matrixVectorMultiply,
    vectorNorm,
    vectorSubtract,
    validateExpression,
    getDiagonal,
    getLowerTriangular,
    getUpperTriangular,
    matrixAdd,
    matrixInverse,
    matrixMultiply,
    eigenvalues,
    spectralRadius
  };
}
