/**
 * MÉTODOS DE INTERPOLACIÓN
 * Interpolación polinomial y splines
 */

// ============================================
// INTERPOLACIÓN POLINOMIAL
// ============================================

/**
 * INTERPOLACIÓN DE VANDERMONDE
 */
function vandermonde(X, Y) {
  const n = X.length;
  
  // Construir matriz de Vandermonde
  const A = [];
  for (let i = 0; i < n; i++) {
    const row = [];
    for (let j = 0; j < n; j++) {
      row.push(Math.pow(X[i], n - 1 - j));
    }
    A.push(row);
  }
  
  // Resolver sistema A * coef = Y para obtener coeficientes
  const coef = partialPivoting(A, Y, false);
  
  return {
    coefficients: coef,
    degree: n - 1
  };
}

/**
 * DIFERENCIAS DIVIDIDAS DE NEWTON
 */
function dividedDifferences(X, Y) {
  const n = X.length;
  
  // Inicializar tabla de diferencias divididas
  const D = Array(n).fill(0).map(() => Array(n).fill(0));
  
  // Primera columna es Y
  for (let i = 0; i < n; i++) {
    D[i][0] = Y[i];
  }
  
  // Calcular diferencias divididas
  for (let j = 1; j < n; j++) {
    for (let i = j; i < n; i++) {
      D[i][j] = (D[i][j - 1] - D[i - 1][j - 1]) / (X[i] - X[i - j]);
    }
  }
  
  // Coeficientes están en la diagonal
  const coef = [];
  for (let i = 0; i < n; i++) {
    coef.push(D[i][i]);
  }
  
  return {
    coefficients: coef,
    x_points: X,
    degree: n - 1,
    table: D
  };
}

/**
 * INTERPOLACIÓN DE LAGRANGE
 */
function lagrange(X, Y) {
  const n = X.length;
  
  // Almacenar polinomios base de Lagrange
  const L = [];
  
  // Calcular cada polinomio base L_i(x)
  for (let i = 0; i < n; i++) {
    // Empezar con polinomio [1]
    let Li = [1.0];
    
    // Multiplicar por (x - X[j]) / (X[i] - X[j]) para todos j != i
    for (let j = 0; j < n; j++) {
      if (i !== j) {
        // Polinomio (x - X[j])
        const numerator = [1.0, -X[j]];
        // Dividir por (X[i] - X[j])
        const denominator = X[i] - X[j];
        // Multiplicar Li por este factor
        Li = convolve(Li, numerator).map(val => val / denominator);
      }
    }
    
    L.push(Li);
  }
  
  // Calcular polinomio final P(x) = sum(Y[i] * L_i(x))
  let P = new Array(n).fill(0);
  for (let i = 0; i < n; i++) {
    P = polyadd(P, L[i].map(val => Y[i] * val));
  }
  
  return {
    coefficients: P,
    lagrange_basis: L,
    x_points: X,
    y_points: Y,
    degree: n - 1
  };
}

// ============================================
// SPLINES (TRAZADORES)
// ============================================

/**
 * SPLINE LINEAL (Trazadores lineales)
 * Cada segmento: S_i(x) = a_i*x + b_i
 */
function linearSpline(X, Y) {
  const n = X.length;
  
  const Coef = [];
  for (let i = 0; i < n - 1; i++) {
    // Verificar división por cero
    if (X[i + 1] === X[i]) {
      throw new Error(`Duplicate x-values at index ${i} and ${i + 1}: X[${i}] = X[${i + 1}] = ${X[i]}. X values must be distinct.`);
    }
    
    // Calcular pendiente (a)
    const a_i = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
    
    // Calcular intersección (b)
    const b_i = Y[i] - a_i * X[i];
    
    Coef.push([a_i, b_i]);
  }
  
  return {
    coefficients: Coef,
    x_points: X,
    y_points: Y,
    num_segments: n - 1
  };
}

/**
 * SPLINE CUADRÁTICO (Trazadores cuadráticos)
 * Cada segmento: S_i(x) = a_i*x^2 + b_i*x + c_i
 * Usa condición: a_0 = 0 (primer segmento es lineal)
 */
function quadraticSpline(X, Y) {
  const n = X.length;
  
  // Verificar mínimo 3 puntos
  if (n < 3) {
    throw new Error('Quadratic spline requires at least 3 points');
  }
  
  // Verificar valores X duplicados y ordenados
  for (let i = 0; i < n - 1; i++) {
    if (X[i + 1] <= X[i]) {
      throw new Error(`X values must be strictly increasing. Problem at X[${i}]=${X[i]}, X[${i+1}]=${X[i+1]}`);
    }
  }
  
  const m = 3 * (n - 1); // Número de incógnitas
  
  // Inicializar sistema de ecuaciones A*coef = b
  const A = Array(m).fill(0).map(() => Array(m).fill(0));
  const b = Array(m).fill(0);
  
  let row = 0;
  
  // Condiciones de interpolación: S_i(x_i) = y_i y S_i(x_{i+1}) = y_{i+1}
  // Total: 2(n-1) ecuaciones
  for (let i = 0; i < n - 1; i++) {
    const colStart = 3 * i;
    
    // S_i(x_i) = y_i
    A[row][colStart] = X[i] * X[i];
    A[row][colStart + 1] = X[i];
    A[row][colStart + 2] = 1;
    b[row] = Y[i];
    row++;
    
    // S_i(x_{i+1}) = y_{i+1}
    A[row][colStart] = X[i + 1] * X[i + 1];
    A[row][colStart + 1] = X[i + 1];
    A[row][colStart + 2] = 1;
    b[row] = Y[i + 1];
    row++;
  }
  
  // Condiciones de continuidad de primera derivada: S'_{i-1}(x_i) = S'_i(x_i)
  // Total: (n-2) ecuaciones
  for (let i = 1; i < n - 1; i++) {
    const colStartPrev = 3 * (i - 1);
    const colStartCurr = 3 * i;
    
    // S'_{i-1}(x_i) = 2*a_{i-1}*x_i + b_{i-1}
    // S'_i(x_i) = 2*a_i*x_i + b_i
    A[row][colStartPrev] = 2 * X[i];     // a_{i-1}
    A[row][colStartPrev + 1] = 1;        // b_{i-1}
    A[row][colStartCurr] = -2 * X[i];    // -a_i
    A[row][colStartCurr + 1] = -1;       // -b_i
    b[row] = 0;
    row++;
  }
  
  // Condición de frontera: a_0 = 0 (primer segmento es lineal)
  // Esta es la ecuación adicional necesaria
  A[row][0] = 1;
  b[row] = 0;
  row++;
  
  // Verificar que tenemos el número correcto de ecuaciones
  if (row !== m) {
    throw new Error(`System mismatch: ${row} equations for ${m} unknowns`);
  }
  
  // Resolver el sistema usando pivoteo parcial para mayor estabilidad
  let Saux;
  try {
    Saux = partialPivoting(A, b);
  } catch (e) {
    // Si falla con pivoteo, intentar Gauss simple
    console.warn('Partial pivoting failed, trying simple Gauss elimination');
    Saux = gaussElimination(A, b);
  }
  
  // Organizar coeficientes: cada segmento tiene [a_i, b_i, c_i]
  const Coef = [];
  for (let i = 0; i < n - 1; i++) {
    Coef.push([
      Saux[3 * i],      // a_i
      Saux[3 * i + 1],  // b_i
      Saux[3 * i + 2]   // c_i
    ]);
  }
  
  return {
    coefficients: Coef,
    x_points: X,
    y_points: Y,
    num_segments: n - 1
  };
}

/**
 * SPLINE CÚBICO (Trazadores cúbicos naturales)
 * Cada segmento: S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
 * Usando forma de Hermite para mejor condicionamiento
 */
function cubicSpline(X, Y) {
  const n = X.length;
  
  // Verificar mínimo 3 puntos
  if (n < 3) {
    throw new Error('Cubic spline requires at least 3 points');
  }
  
  // Verificar valores X ordenados
  for (let i = 0; i < n - 1; i++) {
    if (X[i + 1] <= X[i]) {
      throw new Error(`X values must be strictly increasing. Problem at X[${i}]=${X[i]}, X[${i+1}]=${X[i+1]}`);
    }
  }
  
  // Método simplificado usando segunda derivada
  // Construir sistema tridiagonal para M_i (segunda derivada en nodos)
  const h = []; // h[i] = X[i+1] - X[i]
  for (let i = 0; i < n - 1; i++) {
    h.push(X[i + 1] - X[i]);
  }
  
  // Construir sistema tridiagonal A*M = b para segundas derivadas
  const A = Array(n).fill(0).map(() => Array(n).fill(0));
  const b = Array(n).fill(0);
  
  // Primera ecuación: M_0 = 0 (frontera natural)
  A[0][0] = 1;
  b[0] = 0;
  
  // Ecuaciones intermedias
  for (let i = 1; i < n - 1; i++) {
    A[i][i - 1] = h[i - 1];
    A[i][i] = 2 * (h[i - 1] + h[i]);
    A[i][i + 1] = h[i];
    
    b[i] = 6 * ((Y[i + 1] - Y[i]) / h[i] - (Y[i] - Y[i - 1]) / h[i - 1]);
  }
  
  // Última ecuación: M_{n-1} = 0 (frontera natural)
  A[n - 1][n - 1] = 1;
  b[n - 1] = 0;
  
  // Resolver el sistema para obtener M (segundas derivadas)
  let M;
  try {
    M = partialPivoting(A, b);
  } catch (e) {
    console.warn('Partial pivoting failed, trying simple Gauss elimination');
    M = gaussElimination(A, b);
  }
  
  // Calcular coeficientes de cada segmento
  // S_i(x) = a_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3
  const Coef = [];
  
  for (let i = 0; i < n - 1; i++) {
    const a_i = Y[i];
    const b_i = (Y[i + 1] - Y[i]) / h[i] - h[i] * (2 * M[i] + M[i + 1]) / 6;
    const c_i = M[i] / 2;
    const d_i = (M[i + 1] - M[i]) / (6 * h[i]);
    
    Coef.push([a_i, b_i, c_i, d_i]);
  }
  
  return {
    coefficients: Coef,
    x_points: X,
    y_points: Y,
    num_segments: n - 1,
    second_derivatives: M
  };
}

// ============================================
// FUNCIONES AUXILIARES
// ============================================

/**
 * Convolución de dos polinomios (multiplicación)
 * @param {Array<number>} a - Coeficientes del primer polinomio
 * @param {Array<number>} b - Coeficientes del segundo polinomio
 * @returns {Array<number>} - Coeficientes del producto
 */
function convolve(a, b) {
  const result = new Array(a.length + b.length - 1).fill(0);
  
  for (let i = 0; i < a.length; i++) {
    for (let j = 0; j < b.length; j++) {
      result[i + j] += a[i] * b[j];
    }
  }
  
  return result;
}

/**
 * Suma de dos polinomios
 * @param {Array<number>} a - Coeficientes del primer polinomio
 * @param {Array<number>} b - Coeficientes del segundo polinomio
 * @returns {Array<number>} - Coeficientes de la suma
 */
function polyadd(a, b) {
  const maxLen = Math.max(a.length, b.length);
  const result = new Array(maxLen).fill(0);
  
  for (let i = 0; i < a.length; i++) {
    result[i] += a[i];
  }
  
  for (let i = 0; i < b.length; i++) {
    result[i] += b[i];
  }
  
  return result;
}

// Exportar funciones
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    vandermonde,
    dividedDifferences,
    lagrange,
    linearSpline,
    quadraticSpline,
    cubicSpline
  };
}
