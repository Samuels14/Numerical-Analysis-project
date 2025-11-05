/**
 * MÉTODOS PARA SISTEMAS LINEALES
 * Métodos directos e iterativos para resolver Ax = b
 */

// ============================================
// MÉTODOS DIRECTOS - ELIMINACIÓN GAUSSIANA
// ============================================

/**
 * ELIMINACIÓN GAUSSIANA SIMPLE
 */
function gaussElimination(A, b, returnStages = false) {
  const n = A.length;
  // Crear matriz aumentada
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  
  if (returnStages) {
    stages.push(['Initial augmented matrix', copyMatrix(M)]);
  }
  
  // Eliminación hacia adelante
  for (let k = 0; k < n - 1; k++) {
    // Verificar pivote cero
    if (Math.abs(M[k][k]) < CONFIG.EPSILON) {
      throw new Error(`Zero pivot encountered at position (${k},${k}). Matrix may be singular.`);
    }
    
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j <= n; j++) {
        M[i][j] -= factor * M[k][j];
      }
    }
    
    if (returnStages) {
      stages.push([`After elimination step ${k + 1}`, copyMatrix(M)]);
    }
  }
  
  // Sustitución hacia atrás
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    if (Math.abs(M[i][i]) < CONFIG.EPSILON) {
      throw new Error(`Zero diagonal element at position (${i},${i}). Matrix may be singular.`);
    }
    
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) {
      s -= M[i][j] * x[j];
    }
    x[i] = s / M[i][i];
  }
  
  if (returnStages) {
    return { solution: x, stages: stages };
  }
  return x;
}

/**
 * PIVOTEO PARCIAL
 */
function partialPivoting(A, b, returnStages = false) {
  const n = A.length;
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  
  if (returnStages) {
    stages.push(['Initial augmented matrix', copyMatrix(M)]);
  }
  
  for (let k = 0; k < n - 1; k++) {
    // Encontrar fila con elemento más grande en columna k
    let maxRow = k;
    for (let i = k + 1; i < n; i++) {
      if (Math.abs(M[i][k]) > Math.abs(M[maxRow][k])) {
        maxRow = i;
      }
    }
    
    // Intercambiar filas
    if (maxRow !== k) {
      [M[k], M[maxRow]] = [M[maxRow], M[k]];
      if (returnStages) {
        stages.push([`After row swap (pivot): row ${k} ↔ row ${maxRow}`, copyMatrix(M)]);
      }
    }
    
    // Eliminación
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j <= n; j++) {
        M[i][j] -= factor * M[k][j];
      }
    }
    
    if (returnStages) {
      stages.push([`After elimination step ${k + 1}`, copyMatrix(M)]);
    }
  }
  
  // Sustitución hacia atrás
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) {
      s -= M[i][j] * x[j];
    }
    x[i] = s / M[i][i];
  }
  
  if (returnStages) {
    return { solution: x, stages: stages };
  }
  return x;
}

/**
 * PIVOTEO TOTAL
 */
function totalPivoting(A, b, returnStages = false) {
  const n = A.length;
  const M = A.map((row, i) => [...row, b[i]]);
  const marks = Array.from({ length: n }, (_, i) => i);
  
  const stages = [];
  
  if (returnStages) {
    stages.push(['Initial augmented matrix', copyMatrix(M)]);
  }
  
  for (let k = 0; k < n - 1; k++) {
    // Encontrar elemento máximo en submatriz
    let maxVal = 0, maxRow = k, maxCol = k;
    for (let i = k; i < n; i++) {
      for (let j = k; j < n; j++) {
        if (Math.abs(M[i][j]) > maxVal) {
          maxVal = Math.abs(M[i][j]);
          maxRow = i;
          maxCol = j;
        }
      }
    }
    
    if (maxVal === 0) continue;
    
    // Intercambiar filas
    if (maxRow !== k) {
      [M[k], M[maxRow]] = [M[maxRow], M[k]];
      if (returnStages) {
        stages.push([`After row swap: row ${k} ↔ row ${maxRow}`, copyMatrix(M)]);
      }
    }
    
    // Intercambiar columnas
    if (maxCol !== k) {
      for (let i = 0; i < n; i++) {
        [M[i][k], M[i][maxCol]] = [M[i][maxCol], M[i][k]];
      }
      [marks[k], marks[maxCol]] = [marks[maxCol], marks[k]];
      if (returnStages) {
        stages.push([`After column swap: col ${k} ↔ col ${maxCol}`, copyMatrix(M)]);
      }
    }
    
    // Eliminación
    for (let i = k + 1; i < n; i++) {
      const factor = M[i][k] / M[k][k];
      for (let j = k; j <= n; j++) {
        M[i][j] -= factor * M[k][j];
      }
    }
    
    if (returnStages) {
      stages.push([`After elimination step ${k + 1}`, copyMatrix(M)]);
    }
  }
  
  // Sustitución hacia atrás
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) {
      s -= M[i][j] * x[j];
    }
    x[i] = s / M[i][i];
  }
  
  // Reordenar solución
  const xFinal = new Array(n);
  for (let i = 0; i < n; i++) {
    xFinal[marks[i]] = x[i];
  }
  
  if (returnStages) {
    return { solution: xFinal, stages: stages, marks: marks };
  }
  return xFinal;
}

// ============================================
// SUSTITUCIÓN PROGRESIVA Y REGRESIVA
// ============================================

/**
 * Sustitución progresiva (forward substitution)
 * Para sistemas triangulares inferiores Lx = b
 */
function forwardSubstitution(L, b) {
  const n = L.length;
  const x = new Array(n).fill(0);
  
  for (let i = 0; i < n; i++) {
    let s = b[i];
    for (let j = 0; j < i; j++) {
      s -= L[i][j] * x[j];
    }
    x[i] = s / L[i][i];
  }
  
  return x;
}

/**
 * Sustitución regresiva (backward substitution)
 * Para sistemas triangulares superiores Ux = b
 */
function backwardSubstitution(U, b) {
  const n = U.length;
  const x = new Array(n).fill(0);
  
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let j = i + 1; j < n; j++) {
      s -= U[i][j] * x[j];
    }
    x[i] = s / U[i][i];
  }
  
  return x;
}

/**
 * Sustitución progresiva para matriz aumentada [L|b]
 * Versión standalone para interfaz de usuario
 */
function forwardSubstitutionMethod(M) {
  const n = M.length;
  const x = new Array(n).fill(0);
  const stages = [];
  
  // Guardar matriz inicial
  stages.push(['Initial augmented matrix [L|b]', copyMatrix(M)]);
  
  // Sustitución progresiva
  const iterationDetails = [];
  x[0] = M[0][n] / M[0][0];
  iterationDetails.push(`x[0] = ${M[0][n].toFixed(6)} / ${M[0][0].toFixed(6)} = ${x[0].toFixed(6)}`);
  
  for (let i = 1; i < n; i++) {
    let s = M[i][n];
    for (let j = 0; j < i; j++) {
      s -= M[i][j] * x[j];
    }
    x[i] = s / M[i][i];
    iterationDetails.push(`x[${i}] = (${M[i][n].toFixed(6)} - sum) / ${M[i][i].toFixed(6)} = ${x[i].toFixed(6)}`);
  }
  
  stages.push(['Forward substitution steps', iterationDetails]);
  
  return { solution: x, stages: stages };
}

/**
 * Sustitución regresiva para matriz aumentada [U|b]
 * Versión standalone para interfaz de usuario
 */
function backwardSubstitutionMethod(M) {
  const n = M.length;
  const x = new Array(n).fill(0);
  const stages = [];
  
  // Guardar matriz inicial
  stages.push(['Initial augmented matrix [U|b]', copyMatrix(M)]);
  
  // Sustitución regresiva
  const iterationDetails = [];
  x[n - 1] = M[n - 1][n] / M[n - 1][n - 1];
  iterationDetails.push(`x[${n - 1}] = ${M[n - 1][n].toFixed(6)} / ${M[n - 1][n - 1].toFixed(6)} = ${x[n - 1].toFixed(6)}`);
  
  for (let i = n - 2; i >= 0; i--) {
    let s = M[i][n];
    for (let j = i + 1; j < n; j++) {
      s -= M[i][j] * x[j];
    }
    x[i] = s / M[i][i];
    iterationDetails.push(`x[${i}] = (${M[i][n].toFixed(6)} - sum) / ${M[i][i].toFixed(6)} = ${x[i].toFixed(6)}`);
  }
  
  stages.push(['Backward substitution steps', iterationDetails]);
  
  return { solution: x, stages: stages };
}

// ============================================
// FACTORIZACIÓN LU
// ============================================

/**
 * FACTORIZACIÓN LU SIMPLE (Gaussiana)
 */
function luSimple(A, b) {
  const n = A.length;
  
  // Inicializar L como identidad y U como ceros
  const L = identityMatrix(n);
  const U = zerosMatrix(n, n);
  
  // Crear matriz aumentada
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Factorización
  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      if (M[j][i] !== 0) {
        L[j][i] = M[j][i] / M[i][i];
        for (let k = i; k <= n; k++) {
          M[j][k] -= L[j][i] * M[i][k];
        }
      }
    }
    
    // Copiar fila a U
    for (let k = i; k < n; k++) {
      U[i][k] = M[i][k];
    }
    
    stages.push({
      M: copyMatrix(M),
      L: copyMatrix(L),
      U: copyMatrix(U)
    });
  }
  
  // Copiar última fila a U
  for (let k = 0; k < n; k++) {
    U[n - 1][k] = M[n - 1][k];
  }
  
  // Resolver Lz = b
  const z = forwardSubstitution(L, b);
  
  // Resolver Ux = z
  const x = backwardSubstitution(U, z);
  
  return { solution: x, L: L, U: U, stages: stages };
}

/**
 * FACTORIZACIÓN LU CON PIVOTEO PARCIAL
 */
function luPartialPivot(A, b) {
  const n = A.length;
  
  const L = identityMatrix(n);
  const U = zerosMatrix(n, n);
  const P = identityMatrix(n);
  
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U),
    P: copyMatrix(P)
  });
  
  for (let i = 0; i < n - 1; i++) {
    // Encontrar fila con elemento más grande
    let maxRow = i;
    let maxVal = Math.abs(M[i][i]);
    
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(M[k][i]) > maxVal) {
        maxVal = Math.abs(M[k][i]);
        maxRow = k;
      }
    }
    
    // Intercambiar filas
    if (maxRow !== i) {
      [M[i], M[maxRow]] = [M[maxRow], M[i]];
      [P[i], P[maxRow]] = [P[maxRow], P[i]];
      
      if (i > 0) {
        for (let k = 0; k < i; k++) {
          [L[i][k], L[maxRow][k]] = [L[maxRow][k], L[i][k]];
        }
      }
    }
    
    // Eliminación
    for (let j = i + 1; j < n; j++) {
      if (M[j][i] !== 0) {
        L[j][i] = M[j][i] / M[i][i];
        for (let k = i; k <= n; k++) {
          M[j][k] -= L[j][i] * M[i][k];
        }
      }
    }
    
    // Copiar fila a U
    for (let k = i; k < n; k++) {
      U[i][k] = M[i][k];
    }
    
    stages.push({
      M: copyMatrix(M),
      L: copyMatrix(L),
      U: copyMatrix(U),
      P: copyMatrix(P)
    });
  }
  
  // Copiar última fila a U
  for (let k = 0; k < n; k++) {
    U[n - 1][k] = M[n - 1][k];
  }
  
  // Calcular Pb
  const Pb = matrixVectorMultiply(P, b);
  
  // Resolver Lz = Pb
  const z = forwardSubstitution(L, Pb);
  
  // Resolver Ux = z
  const x = backwardSubstitution(U, z);
  
  return { solution: x, L: L, U: U, P: P, stages: stages };
}

/**
 * FACTORIZACIÓN DE CROUT
 */
function crout(A, b) {
  const n = A.length;
  
  const L = identityMatrix(n);
  const U = identityMatrix(n);
  
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Factorización
  for (let i = 0; i < n - 1; i++) {
    // Calcular columna i de L (desde fila i hasta n)
    for (let j = i; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[j][k] * U[k][i];
      }
      L[j][i] = A[j][i] - sumVal;
    }
    
    // Calcular fila i de U (desde columna i+1 hasta n)
    for (let j = i + 1; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[i][k] * U[k][j];
      }
      U[i][j] = (A[i][j] - sumVal) / L[i][i];
    }
    
    stages.push({
      M: copyMatrix(M),
      L: copyMatrix(L),
      U: copyMatrix(U)
    });
  }
  
  // Calcular último elemento L[n-1][n-1]
  let sumVal = 0;
  for (let k = 0; k < n - 1; k++) {
    sumVal += L[n - 1][k] * U[k][n - 1];
  }
  L[n - 1][n - 1] = A[n - 1][n - 1] - sumVal;
  
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Resolver Lz = b
  const z = forwardSubstitution(L, b);
  
  // Resolver Ux = z
  const x = backwardSubstitution(U, z);
  
  return { solution: x, L: L, U: U, stages: stages };
}

/**
 * FACTORIZACIÓN DE DOOLITTLE
 */
function doolittle(A, b) {
  const n = A.length;
  
  const L = identityMatrix(n);
  const U = identityMatrix(n);
  
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Factorización
  for (let i = 0; i < n - 1; i++) {
    // Calcular fila i de U (desde columna i hasta n)
    for (let j = i; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[i][k] * U[k][j];
      }
      U[i][j] = A[i][j] - sumVal;
    }
    
    // Calcular columna i de L (desde fila i+1 hasta n)
    for (let j = i + 1; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[j][k] * U[k][i];
      }
      L[j][i] = (A[j][i] - sumVal) / U[i][i];
    }
    
    stages.push({
      M: copyMatrix(M),
      L: copyMatrix(L),
      U: copyMatrix(U)
    });
  }
  
  // Calcular último elemento U[n-1][n-1]
  let sumVal = 0;
  for (let k = 0; k < n - 1; k++) {
    sumVal += L[n - 1][k] * U[k][n - 1];
  }
  U[n - 1][n - 1] = A[n - 1][n - 1] - sumVal;
  
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Resolver Lz = b
  const z = forwardSubstitution(L, b);
  
  // Resolver Ux = z
  const x = backwardSubstitution(U, z);
  
  return { solution: x, L: L, U: U, stages: stages };
}

/**
 * FACTORIZACIÓN DE CHOLESKY
 * Para matrices simétricas definidas positivas
 */
function cholesky(A, b) {
  const n = A.length;
  
  const L = identityMatrix(n);
  const U = identityMatrix(n);
  
  const M = A.map((row, i) => [...row, b[i]]);
  
  const stages = [];
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Factorización
  for (let i = 0; i < n - 1; i++) {
    // Calcular L[i][i] y U[i][i]
    let sumVal = 0;
    for (let k = 0; k < i; k++) {
      sumVal += L[i][k] * U[k][i];
    }
    const value = A[i][i] - sumVal;
    
    if (value < 0) {
      throw new Error(`Matrix is not positive definite (negative value ${value} at position (${i},${i}))`);
    }
    
    L[i][i] = Math.sqrt(value);
    U[i][i] = L[i][i];
    
    // Calcular columna i de L (debajo de la diagonal)
    for (let j = i + 1; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[j][k] * U[k][i];
      }
      L[j][i] = (A[j][i] - sumVal) / U[i][i];
    }
    
    // Calcular fila i de U (a la derecha de la diagonal)
    for (let j = i + 1; j < n; j++) {
      let sumVal = 0;
      for (let k = 0; k < i; k++) {
        sumVal += L[i][k] * U[k][j];
      }
      U[i][j] = (A[i][j] - sumVal) / L[i][i];
    }
    
    stages.push({
      M: copyMatrix(M),
      L: copyMatrix(L),
      U: copyMatrix(U)
    });
  }
  
  // Calcular último elemento diagonal
  let sumVal = 0;
  for (let k = 0; k < n - 1; k++) {
    sumVal += L[n - 1][k] * U[k][n - 1];
  }
  const value = A[n - 1][n - 1] - sumVal;
  
  if (value < 0) {
    throw new Error(`Matrix is not positive definite (negative value ${value} at position (${n - 1},${n - 1}))`);
  }
  
  L[n - 1][n - 1] = Math.sqrt(value);
  U[n - 1][n - 1] = L[n - 1][n - 1];
  
  stages.push({
    M: copyMatrix(M),
    L: copyMatrix(L),
    U: copyMatrix(U)
  });
  
  // Resolver Lz = b
  const z = forwardSubstitution(L, b);
  
  // Resolver Ux = z
  const x = backwardSubstitution(U, z);
  
  return { solution: x, L: L, U: U, stages: stages };
}

// ============================================
// MÉTODOS ITERATIVOS
// ============================================

/**
 * MÉTODO DE JACOBI
 */
function jacobi(A, b, x0, tol, nmax) {
  const n = A.length;
  
  // Verificar elementos diagonales no cero
  for (let i = 0; i < n; i++) {
    if (Math.abs(A[i][i]) < CONFIG.EPSILON) {
      throw new Error(`Diagonal element A[${i}][${i}] is zero. Jacobi method cannot proceed.`);
    }
  }
  
  // Inicialización
  const D = getDiagonal(A);
  const L = getLowerTriangular(A);
  const U = getUpperTriangular(A);
  const LplusU = matrixAdd(L, U);
  
  // Calcular T = D^-1 * (L + U)
  const Dinv = matrixInverse(D);
  const T = matrixMultiply(Dinv, LplusU);
  
  // Calcular C = D^-1 * b
  const C = matrixVectorMultiply(Dinv, b);
  
  // Radio espectral
  const spectralRad = spectralRadius(T);
  
  let xant = [...x0];
  let E = 1000;
  let cont = 0;
  
  const iterations = [];
  
  // Ciclo
  while (E > tol && cont < nmax) {
    const xact = matrixVectorMultiply(T, xant).map((val, i) => val + C[i]);
    E = vectorNorm(vectorSubtract(xant, xact));
    xant = xact;
    cont = cont + 1;
    
    iterations.push([cont - 1, E, ...xact]);
  }
  
  return {
    solution: xant,
    iterations: iterations,
    converged: E < tol,
    T: T,
    C: C,
    spectral_radius: spectralRad
  };
}

/**
 * MÉTODO DE GAUSS-SEIDEL
 */
function gaussSeidel(A, b, x0, tol, nmax) {
  const n = A.length;
  
  // Verificar elementos diagonales no cero
  for (let i = 0; i < n; i++) {
    if (Math.abs(A[i][i]) < CONFIG.EPSILON) {
      throw new Error(`Diagonal element A[${i}][${i}] is zero. Gauss-Seidel method cannot proceed.`);
    }
  }
  
  // Inicialización
  const D = getDiagonal(A);
  const L = getLowerTriangular(A);
  const U = getUpperTriangular(A);
  
  // Calcular T = (D - L)^-1 * U
  const DminusL = matrixAdd(D, L.map(row => row.map(val => -val)));
  const DminusLinv = matrixInverse(DminusL);
  const T = matrixMultiply(DminusLinv, U);
  
  // Calcular C = (D - L)^-1 * b
  const C = matrixVectorMultiply(DminusLinv, b);
  
  // Radio espectral
  const spectralRad = spectralRadius(T);
  
  let xant = [...x0];
  let E = 1000;
  let cont = 0;
  
  const iterations = [];
  
  // Ciclo
  while (E > tol && cont < nmax) {
    const xact = matrixVectorMultiply(T, xant).map((val, i) => val + C[i]);
    E = vectorNorm(vectorSubtract(xant, xact));
    xant = xact;
    cont = cont + 1;
    
    iterations.push([cont - 1, E, ...xact]);
  }
  
  return {
    solution: xant,
    iterations: iterations,
    converged: E < tol,
    T: T,
    C: C,
    spectral_radius: spectralRad
  };
}

/**
 * MÉTODO SOR (Successive Over-Relaxation)
 */
function sor(A, b, x0, omega, tol, nmax) {
  const n = A.length;
  
  // Verificar elementos diagonales no cero
  for (let i = 0; i < n; i++) {
    if (Math.abs(A[i][i]) < CONFIG.EPSILON) {
      throw new Error(`Diagonal element A[${i}][${i}] is zero. SOR method cannot proceed.`);
    }
  }
  
  // Verificar omega en rango válido
  if (omega <= 0 || omega >= 2) {
    throw new Error(`Relaxation factor omega=${omega} is out of range. Must be 0 < omega < 2.`);
  }
  
  // Inicialización
  const D = getDiagonal(A);
  const L = getLowerTriangular(A);
  const U = getUpperTriangular(A);
  
  // Calcular T = (D - ω*L)^-1 * ((1-ω)*D + ω*U)
  const omegaL = L.map(row => row.map(val => omega * val));
  const DminusOmegaL = matrixAdd(D, omegaL.map(row => row.map(val => -val)));
  const DminusOmegaLinv = matrixInverse(DminusOmegaL);
  
  const oneMinusOmegaD = D.map(row => row.map(val => (1 - omega) * val));
  const omegaU = U.map(row => row.map(val => omega * val));
  const rightSide = matrixAdd(oneMinusOmegaD, omegaU);
  
  const T = matrixMultiply(DminusOmegaLinv, rightSide);
  
  // Calcular C = ω * (D - ω*L)^-1 * b
  const Cb = matrixVectorMultiply(DminusOmegaLinv, b);
  const C = Cb.map(val => omega * val);
  
  // Radio espectral
  const spectralRad = spectralRadius(T);
  
  let xant = [...x0];
  let E = 1000;
  let cont = 0;
  
  const iterations = [];
  
  // Ciclo
  while (E > tol && cont < nmax) {
    const xact = matrixVectorMultiply(T, xant).map((val, i) => val + C[i]);
    E = vectorNorm(vectorSubtract(xant, xact));
    xant = xact;
    cont = cont + 1;
    
    iterations.push([cont - 1, E, ...xact]);
  }
  
  return {
    solution: xant,
    iterations: iterations,
    converged: E < tol,
    T: T,
    C: C,
    spectral_radius: spectralRad,
    omega: omega
  };
}

// Exportar funciones
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    gaussElimination,
    partialPivoting,
    totalPivoting,
    forwardSubstitution,
    backwardSubstitution,
    forwardSubstitutionMethod,
    backwardSubstitutionMethod,
    luSimple,
    luPartialPivot,
    crout,
    doolittle,
    cholesky,
    jacobi,
    gaussSeidel,
    sor
  };
}
