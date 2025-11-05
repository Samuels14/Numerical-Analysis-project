/**
 * FORMATEADORES DE RESULTADOS
 * Convierte los resultados numéricos a HTML formateado
 */

/**
 * Formatea resultados de métodos de búsqueda de raíces
 */
function formatRootResult(result, method) {
  let html = '';
  
  switch (method) {
    case 'incremental_search':
      html = formatIncrementalSearch(result);
      break;
      
    case 'bisection':
      html = formatBisection(result);
      break;
      
    case 'false_position':
      html = formatFalsePosition(result);
      break;
      
    case 'fixed_point':
      html = formatFixedPoint(result);
      break;
      
    case 'newton':
      html = formatNewton(result);
      break;
      
    case 'secant':
      html = formatSecant(result);
      break;
      
    case 'newton_multiple':
      html = formatNewtonMultiple(result);
      break;
  }
  
  return html;
}

function formatIncrementalSearch(result) {
  let html = 'Búsqueda incrementales\n\nResultados:\n\n';
  
  if (result.intervals.length > 0) {
    for (const interval of result.intervals) {
      const [a, b] = interval;
      html += `Hay una raiz de f en [${a.toFixed(10)},${b.toFixed(10)}]\n`;
    }
    html += '_______________________________________________________';
  } else {
    html = 'No se encontraron raíces en el intervalo especificado';
  }
  
  return html;
}

function formatBisection(result) {
  let html = 'Bisección\n\nTabla de resultados:\n\n';
  html += '| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\n';
  
  let prevXm = null;
  for (const iterData of result.iterations) {
    const [i, a, b, c, fc] = iterData;
    let errorStr = '';
    
    if (prevXm !== null) {
      const error = Math.abs(c - prevXm);
      errorStr = `|  ${error.toExponential(1)}  |`;
    } else {
      errorStr = '|           |';
    }
    
    html += `|  ${(i + 1).toString().padEnd(2)} | ${a.toFixed(10)}  | ${c.toFixed(10)}  |  ${b.toFixed(10)}  |  ${fc.toExponential(1)}  ${errorStr} \n`;
    prevXm = c;
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

function formatFalsePosition(result) {
  let html = 'Regla falsa\n\nTabla de resultados:\n\n';
  html += '| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\n';
  
  let prevXm = null;
  for (const iterData of result.iterations) {
    const [i, a, b, c, fc] = iterData;
    let errorStr = '';
    
    if (prevXm !== null) {
      const error = Math.abs(c - prevXm);
      errorStr = `|  ${error.toExponential(1)}  |`;
    } else {
      errorStr = '|           |';
    }
    
    html += `|  ${(i + 1).toString().padEnd(2)} | ${a.toFixed(10)}  | ${c.toFixed(10)}  |  ${b.toFixed(10)}  |  ${fc.toExponential(1)}  ${errorStr} \n`;
    prevXm = c;
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

function formatFixedPoint(result) {
  let html = 'Punto fijo\n\nTabla de resultados:\n\n';
  html += '| iter|     xi         |     g(xi)      |   f(xi)   |     E     |\n';
  
  for (const iterData of result.iterations) {
    const [i, xOld, xNew] = iterData;
    const fx = result.f(xNew);
    const error = Math.abs(xNew - xOld);
    const errorStr = i > 0 ? `|  ${error.toExponential(1)}  |` : '|           |';
    
    html += `|  ${i.toString().padEnd(2)} | ${xOld.toFixed(10)}  | ${xNew.toFixed(10)}  |  ${fx.toExponential(1)}  ${errorStr} \n`;
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

function formatNewton(result) {
  const x0 = result.iterations.length > 0 ? result.iterations[0][1] : 0;
  const fx0 = result.f(x0);
  
  let html = 'Newton\n\nTabla de resultados:\n\n';
  html += '| iter|     xi         |   f(xi)   |     E     |\n';
  html += `|  0  | ${x0.toFixed(10)}   |  ${fx0.toExponential(1)}  |           |  \n`;
  
  for (const iterData of result.iterations) {
    const [i, xVal, fxVal, dfxVal, xNext] = iterData;
    const error = Math.abs(xNext - xVal);
    html += `|  ${(i + 1).toString().padEnd(2)} | ${xNext.toFixed(10)}  |  ${result.f(xNext).toExponential(1)}  |  ${error.toExponential(1)}  | \n`;
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

function formatSecant(result) {
  let html = 'Secante\n\nTabla de resultados:\n\n';
  html += '| iter|     xi         |   f(xi)   |     E     |\n';
  
  if (result.iterations.length > 0) {
    const x0 = result.iterations[0][1];
    const x1 = result.iterations[0][2];
    html += `|  0  | ${x0.toFixed(10)}   |  ${result.f(x0).toExponential(1)}  |           |  \n`;
    html += `|  1  | ${x1.toFixed(10)}   |  ${result.f(x1).toExponential(1)}  |           |  \n`;
    
    let prevX = x1;
    for (let idx = 0; idx < result.iterations.length; idx++) {
      const iterData = result.iterations[idx];
      const [i, x0Val, x1Val, x2Val, fx1Val] = iterData;
      const error = Math.abs(x2Val - prevX);
      html += `|  ${(idx + 2).toString().padEnd(2)} | ${x2Val.toFixed(10)}  |  ${result.f(x2Val).toExponential(1)}  |  ${error.toExponential(1)}  | \n`;
      prevX = x2Val;
    }
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

function formatNewtonMultiple(result) {
  const x0 = result.iterations.length > 0 ? result.iterations[0][1] : 0;
  const fx0 = result.f(x0);
  
  let html = 'Raíces Múltiples\n\nTabla de resultados:\n\n';
  html += '| iter|     xi         |   f(xi)   |     E     |\n';
  html += `|  0  | ${x0.toFixed(10)}   |  ${fx0.toExponential(1)}  |           |  \n`;
  
  for (const iterData of result.iterations) {
    const [i, xVal, fxVal, dfxVal, d2fxVal, xNext] = iterData;
    const error = Math.abs(xNext - xVal);
    html += `|  ${(i + 1).toString().padEnd(2)} | ${xNext.toFixed(10)}   |  ${result.f(xNext).toExponential(1)}  |  ${error.toExponential(1)}  | \n`;
  }
  
  html += `\nSe encontró una aproximación de la raiz en ${result.root.toFixed(15)}`;
  
  return html;
}

/**
 * Formatea resultados de sistemas lineales
 */
function formatLinearResult(result, method) {
  let html = '';
  
  // Métodos con etapas
  if (result.stages && Array.isArray(result.stages)) {
    if (method === 'forward_substitution' || method === 'backward_substitution') {
      html = formatSubstitutionResult(result, method);
    } else if (['gauss', 'partial_pivot', 'total_pivot'].includes(method)) {
      html = formatEliminationResult(result, method);
    } else {
      html = formatFactorizationResult(result, method);
    }
  }
  
  // Métodos iterativos
  else if (result.iterations) {
    html = formatIterativeResult(result, method);
  }
  
  return html;
}

function formatEliminationResult(result, method) {
  const methodNames = {
    'gauss': 'Eliminación gaussiana simple',
    'partial_pivot': 'Eliminación gaussiana con pivoteo parcial',
    'total_pivot': 'Eliminación gaussiana con pivoteo total'
  };
  
  let html = `${methodNames[method]}\n\nResultados:\n`;
  
  for (let i = 0; i < result.stages.length; i++) {
    const [desc, matrix] = result.stages[i];
    html += `\nEtapa ${i}\n\n`;
    html += formatMatrix(matrix) + '\n';
  }
  
  html += '\n\nDespués de aplicar sustitución regresiva\n\nx:\n';
  for (const x of result.solution) {
    html += `${x.toFixed(6)}\n`;
  }
  
  return html;
}

function formatSubstitutionResult(result, method) {
  const methodName = method === 'forward_substitution' 
    ? 'Sustitución Progresiva' 
    : 'Sustitución Regresiva';
    
  let html = `${methodName}\n\nResultados:\n\n`;
  
  // Mostrar matriz inicial
  const [desc, matrix] = result.stages[0];
  html += `${desc}:\n\n`;
  html += formatMatrix(matrix) + '\n';
  
  // Mostrar pasos de cálculo
  html += '\nPasos de cálculo:\n\n';
  const [desc2, steps] = result.stages[1];
  for (const step of steps) {
    html += step + '\n';
  }
  
  html += '\nSolución:\n\nx:\n';
  for (const xVal of result.solution) {
    html += `${xVal.toFixed(6)}\n`;
  }
  
  return html;
}

function formatFactorizationResult(result, method) {
  const methodNames = {
    'lu_simple': 'LU con gaussiana simple:',
    'lu_partial': 'LU con pivoteo parcial',
    'crout': 'Crout',
    'doolittle': 'Doolittle',
    'cholesky': 'Cholesky'
  };
  
  let html = `${methodNames[method]}\n\nResultados:\n`;
  
  for (let idx = 0; idx < result.stages.length; idx++) {
    const stage = result.stages[idx];
    html += `\nEtapa ${idx}\n\n`;
    html += formatMatrix(stage.M) + '\n';
    
    if (idx > 0) {
      html += '\nL:\n';
      html += formatMatrix(stage.L) + '\n';
      html += '\nU:\n';
      html += formatMatrix(stage.U) + '\n';
      
      if (stage.P) {
        html += '\nP:\n';
        html += formatMatrix(stage.P) + '\n';
      }
    }
  }
  
  html += '\n\nDespués de aplicar sustitución progresiva y regresiva\n\nx:\n';
  for (const x of result.solution) {
    html += `${x.toFixed(6)}\n`;
  }
  
  return html;
}

function formatIterativeResult(result, method) {
  const methodNames = {
    'jacobi': 'Jacobi',
    'gauss_seidel': 'Gauss-Seidel',
    'sor': 'SOR (relajación)'
  };
  
  let html = `${methodNames[method]}\n\nResultados:\n\nT:\n`;
  html += formatMatrix(result.T) + '\n';
  html += '\nC:\n';
  for (const cVal of result.C) {
    html += ` ${cVal.toFixed(6)} `;
  }
  html += '\n\nradio espectral:\n';
  html += ` ${result.spectral_radius.toFixed(6)}\n\n`;
  
  // Tabla de iteraciones
  html += '| iter|      E     |';
  for (let i = 0; i < result.solution.length; i++) {
    html += ' ';
  }
  html += '\n';
  
  // Primera iteración (x0)
  const firstIter = result.iterations[0];
  const x0 = firstIter.slice(2);
  html += `|  0  |            | ${x0.map(x => x.toFixed(6)).join('  ')} \n`;
  
  // Resto de iteraciones
  for (const iterData of result.iterations) {
    const k = iterData[0];
    const error = iterData[1];
    const xVals = iterData.slice(2);
    html += `|  ${(k + 1).toString().padEnd(2)} |   ${error.toExponential(1)}  | ${xVals.map(x => x.toFixed(6)).join('  ')} \n`;
  }
  
  return html;
}

/**
 * Formatea resultados de interpolación
 */
function formatInterpolationResult(result, method) {
  let html = '';
  
  switch (method) {
    case 'vandermonde':
      html = formatVandermonde(result);
      break;
      
    case 'divided_differences':
      html = formatDividedDifferences(result);
      break;
      
    case 'lagrange':
      html = formatLagrange(result);
      break;
      
    case 'linear_spline':
      html = formatLinearSpline(result);
      break;
      
    case 'quadratic_spline':
      html = formatQuadraticSpline(result);
      break;
      
    case 'cubic_spline':
      html = formatCubicSpline(result);
      break;
  }
  
  return html;
}

function formatVandermonde(result) {
  const n = result.coefficients.length;
  const X = result.X || [];
  
  // Construir matriz de Vandermonde
  const V = [];
  for (let i = 0; i < n; i++) {
    const row = [];
    for (let j = 0; j < n; j++) {
      row.push(Math.pow(X[i], n - 1 - j));
    }
    V.push(row);
  }
  
  let html = 'Vandermonde\n\nResultados:\n\nMatriz de Vandermonde:\n\n';
  html += formatMatrix(V) + '\n';
  html += '\nCoeficientes del polinomio:\n\n';
  for (const c of result.coefficients) {
    html += ` ${c.toFixed(6)} `;
  }
  html += '\n\nPolinomio:\n\n';
  
  // Formatear polinomio
  const polyTerms = [];
  for (let i = 0; i < result.coefficients.length; i++) {
    const exp = result.degree - i;
    if (exp > 1) {
      polyTerms.push(`${result.coefficients[i].toFixed(6)}x^${exp}`);
    } else if (exp === 1) {
      polyTerms.push(`${result.coefficients[i].toFixed(6)}x`);
    } else {
      polyTerms.push(`${result.coefficients[i].toFixed(6)}`);
    }
  }
  
  html += polyTerms.join('+').replace(/\+-/g, '-');
  
  return html;
}

function formatDividedDifferences(result) {
  let html = 'Newton\n\nResultados:\n\nTabla de diferencias divididas:\n\n';
  html += formatMatrix(result.table) + '\n';
  html += '\nCoeficientes del polinomio de Newton:\n\n';
  for (const c of result.coefficients) {
    html += ` ${c.toFixed(6)} `;
  }
  html += '\n\nPolinomio de Newton:\n\n';
  
  // Formatear polinomio de Newton
  const polyParts = [`${result.coefficients[0].toFixed(6)}`];
  for (let i = 1; i < result.coefficients.length; i++) {
    let factors = '';
    for (let j = 0; j < i; j++) {
      const sign = result.x_points[j] >= 0 ? '-' : '+';
      const val = Math.abs(result.x_points[j]);
      factors += `(x${sign}${val.toFixed(6)})`;
    }
    polyParts.push(`${result.coefficients[i].toFixed(6)}${factors}`);
  }
  
  html += polyParts.join('+').replace(/\+-/g, '-');
  
  return html;
}

function formatLagrange(result) {
  let html = 'Lagrange\n\nResultados:\n\nPolinomios interpolantes de Lagrange:\n\n';
  
  // Formatear polinomios base de Lagrange
  for (let i = 0; i < result.lagrange_basis.length; i++) {
    const Li = result.lagrange_basis[i];
    const liDegree = Li.length - 1;
    const terms = [];
    
    for (let j = 0; j < Li.length; j++) {
      const exp = liDegree - j;
      if (Math.abs(Li[j]) > 1e-10) {
        if (exp > 1) {
          terms.push(`${Li[j].toFixed(6)}x^${exp}`);
        } else if (exp === 1) {
          terms.push(`${Li[j].toFixed(6)}x`);
        } else {
          terms.push(`${Li[j].toFixed(6)}`);
        }
      }
    }
    html += terms.join('+').replace(/\+-/g, '-') + `   //L${i}\n`;
  }
  
  html += '\n\nPolinomio:\n\n';
  const polyParts = [];
  for (let i = 0; i < result.y_points.length; i++) {
    polyParts.push(`${result.y_points[i].toFixed(1)}*L${i}`);
  }
  html += polyParts.join('+').replace(/\+-/g, '-');
  
  return html;
}

function formatLinearSpline(result) {
  let html = 'Trazadores lineales\n\nResultados:\n\nCoeficientes de los trazadores:\n\n';
  
  for (const [a, b] of result.coefficients) {
    html += ` ${a.toFixed(6)}  ${b.toFixed(6)} \n`;
  }
  
  html += '\nTrazadores:\n\n';
  for (const [a, b] of result.coefficients) {
    const sign = b >= 0 ? '+' : '';
    html += `${a.toFixed(6)}x${sign}${b.toFixed(6)}\n`;
  }
  
  return html;
}

function formatQuadraticSpline(result) {
  let html = 'Trazadores cuadráticos\n\nResultados:\n\nCoeficientes de los trazadores:\n\n';
  
  for (const [a, b, c] of result.coefficients) {
    html += ` ${a.toFixed(6)}  ${b.toFixed(6)}  ${c.toFixed(6)} \n`;
  }
  
  html += '\nTrazadores:\n\n';
  for (const [a, b, c] of result.coefficients) {
    const bSign = b >= 0 ? '+' : '';
    const cSign = c >= 0 ? '+' : '';
    html += `${a.toFixed(6)}x^2${bSign}${b.toFixed(6)}x${cSign}${c.toFixed(6)}\n`;
  }
  
  return html;
}

function formatCubicSpline(result) {
  let html = '<h4>Spline Cúbico Natural</h4>';
  
  html += '<p><strong>Forma de cada segmento:</strong> S<sub>i</sub>(x) = a<sub>i</sub> + b<sub>i</sub>·(x-x<sub>i</sub>) + c<sub>i</sub>·(x-x<sub>i</sub>)² + d<sub>i</sub>·(x-x<sub>i</sub>)³</p>';
  
  html += '<h5>Coeficientes:</h5>';
  html += '<table>';
  html += '<tr><th>Segmento</th><th>Intervalo</th><th>a<sub>i</sub></th><th>b<sub>i</sub></th><th>c<sub>i</sub></th><th>d<sub>i</sub></th></tr>';
  
  for (let i = 0; i < result.coefficients.length; i++) {
    const [a, b, c, d] = result.coefficients[i];
    const x0 = result.x_points[i];
    const x1 = result.x_points[i + 1];
    
    html += '<tr>';
    html += `<td>S<sub>${i}</sub>(x)</td>`;
    html += `<td>[${x0.toFixed(2)}, ${x1.toFixed(2)}]</td>`;
    html += `<td>${a.toFixed(6)}</td>`;
    html += `<td>${b.toFixed(6)}</td>`;
    html += `<td>${c.toFixed(6)}</td>`;
    html += `<td>${d.toFixed(6)}</td>`;
    html += '</tr>';
  }
  html += '</table>';
  
  html += '<h5>Polinomios por segmento:</h5>';
  for (let i = 0; i < result.coefficients.length; i++) {
    const [a, b, c, d] = result.coefficients[i];
    const x0 = result.x_points[i];
    const x1 = result.x_points[i + 1];
    
    const bSign = b >= 0 ? '+' : '';
    const cSign = c >= 0 ? '+' : '';
    const dSign = d >= 0 ? '+' : '';
    
    html += `<p><strong>S<sub>${i}</sub>(x)</strong> = ${a.toFixed(4)} ${bSign} ${b.toFixed(4)}·(x-${x0}) ${cSign} ${c.toFixed(4)}·(x-${x0})² ${dSign} ${d.toFixed(4)}·(x-${x0})³</p>`;
    html += `<p style="margin-left: 20px; color: #666;">para x ∈ [${x0.toFixed(2)}, ${x1.toFixed(2)}]</p>`;
  }
  
  return html;
}
