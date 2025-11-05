/**
 * MÉTODOS PARA BÚSQUEDA DE RAÍCES DE ECUACIONES
 * Equivalente a los métodos de Python para resolver f(x) = 0
 */

/**
 * BÚSQUEDA INCREMENTAL
 * Encuentra todos los intervalos donde hay cambio de signo
 */
function incrementalSearch(f, x0, delta, nmax) {
  const iterations = [];
  const intervals = [];
  
  // Inicialización
  let xant = x0;
  let fant = f(xant);
  let xact = xant + delta;
  let fact = f(xact);
  
  // Ciclo - buscar TODOS los intervalos
  for (let i = 0; i < nmax; i++) {
    iterations.push([i, xant, fant, xact, fact]);
    
    // Verificar cambio de signo
    if (fant * fact < 0) {
      intervals.push([xant, xact]);
    }
    
    // Moverse al siguiente intervalo
    xant = xact;
    fant = fact;
    xact = xant + delta;
    fact = f(xact);
  }
  
  return {
    intervals: intervals,
    iterations: iterations
  };
}

/**
 * MÉTODO DE BISECCIÓN
 * Criterio de parada: E = |pm - p0| < tol
 */
function bisectionMethod(f, a, b, tol, nmax) {
  const iterations = [];
  
  // Inicialización
  let fa = f(a);
  let pm = (a + b) / 2;
  let fpm = f(pm);
  let E = 1000;
  let cont = 1;
  
  iterations.push([0, a, b, pm, fpm]);
  
  // Ciclo
  while (E > tol && cont < nmax) {
    if (fa * fpm < 0) {
      b = pm;
    } else {
      a = pm;
    }
    
    const p0 = pm;
    pm = (a + b) / 2;
    fpm = f(pm);
    E = Math.abs(pm - p0);
    
    iterations.push([cont, a, b, pm, fpm]);
    cont = cont + 1;
    fa = f(a);
  }
  
  return {
    root: pm,
    iterations: iterations,
    error: E,
    iter_count: cont
  };
}

/**
 * MÉTODO DE REGLA FALSA (FALSE POSITION)
 * Criterio de parada: E = |pm - p0| < tol
 */
function falsePositionMethod(f, a, b, tol, nmax) {
  const iterations = [];
  
  // Inicialización
  let fa = f(a);
  let fb = f(b);
  let pm = (fb * a - fa * b) / (fb - fa);
  let fpm = f(pm);
  let E = 1000;
  let cont = 1;
  
  iterations.push([0, a, b, pm, fpm]);
  
  // Ciclo
  while (E > tol && cont < nmax) {
    if (fa * fpm < 0) {
      b = pm;
    } else {
      a = pm;
    }
    
    const p0 = pm;
    fa = f(a);
    fb = f(b);
    pm = (fb * a - fa * b) / (fb - fa);
    fpm = f(pm);
    E = Math.abs(pm - p0);
    
    iterations.push([cont, a, b, pm, fpm]);
    cont = cont + 1;
  }
  
  return {
    root: pm,
    iterations: iterations,
    error: E,
    iter_count: cont
  };
}

/**
 * MÉTODO DE PUNTO FIJO
 */
function fixedPointMethod(g, x0, tol, nmax) {
  const iterations = [];
  
  for (let i = 0; i < nmax; i++) {
    const x1 = g(x0);
    iterations.push([i, x0, x1]);
    
    if (Math.abs(x1 - x0) < tol) {
      return {
        root: x1,
        iterations: iterations
      };
    }
    
    x0 = x1;
  }
  
  return {
    root: x0,
    iterations: iterations
  };
}

/**
 * MÉTODO DE NEWTON-RAPHSON
 */
function newtonMethod(f, df, x0, tol, nmax) {
  const iterations = [];
  
  for (let i = 0; i < nmax; i++) {
    const fx = f(x0);
    const dfx = df ? df(x0) : numericDerivative(f, x0);
    
    if (Math.abs(dfx) < CONFIG.EPSILON) {
      break;
    }
    
    const x1 = x0 - fx / dfx;
    iterations.push([i, x0, fx, dfx, x1]);
    
    if (Math.abs(x1 - x0) < tol) {
      return {
        root: x1,
        iterations: iterations
      };
    }
    
    x0 = x1;
  }
  
  return {
    root: x0,
    iterations: iterations
  };
}

/**
 * MÉTODO DE LA SECANTE
 */
function secantMethod(f, x0, x1, tol, nmax) {
  const iterations = [];
  
  for (let i = 0; i < nmax; i++) {
    const fx0 = f(x0);
    const fx1 = f(x1);
    
    if (Math.abs(fx1 - fx0) < CONFIG.EPSILON) {
      break;
    }
    
    const x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0);
    iterations.push([i, x0, x1, x2, fx1]);
    
    if (Math.abs(x2 - x1) < tol) {
      return {
        root: x2,
        iterations: iterations
      };
    }
    
    x0 = x1;
    x1 = x2;
  }
  
  return {
    root: x1,
    iterations: iterations
  };
}

/**
 * MÉTODO DE NEWTON PARA RAÍCES MÚLTIPLES
 */
function newtonMultipleRoots(f, df, d2f, x0, tol, nmax) {
  const iterations = [];
  
  for (let i = 0; i < nmax; i++) {
    const fx = f(x0);
    const dfx = df ? df(x0) : numericDerivative(f, x0);
    const d2fx = d2f ? d2f(x0) : numericSecondDerivative(f, x0);
    
    const denom = dfx ** 2 - fx * d2fx;
    if (Math.abs(denom) < CONFIG.EPSILON) {
      break;
    }
    
    const x1 = x0 - (fx * dfx) / denom;
    iterations.push([i, x0, fx, dfx, d2fx, x1]);
    
    if (Math.abs(x1 - x0) < tol) {
      return {
        root: x1,
        iterations: iterations
      };
    }
    
    x0 = x1;
  }
  
  return {
    root: x0,
    iterations: iterations
  };
}

// Exportar funciones
if (typeof module !== 'undefined' && module.exports) {
  module.exports = {
    incrementalSearch,
    bisectionMethod,
    falsePositionMethod,
    fixedPointMethod,
    newtonMethod,
    secantMethod,
    newtonMultipleRoots
  };
}
