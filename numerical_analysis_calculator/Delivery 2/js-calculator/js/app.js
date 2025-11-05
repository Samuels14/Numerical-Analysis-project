/**
 * CONTROLADOR PRINCIPAL
 * Conecta la interfaz con los métodos numéricos
 */

// Estado global de la aplicación
const AppState = {
  currentMethod: 'bisection',
  lastResult: null,
  calculator: null // Para Desmos/GeoGebra
};

/**
 * Inicializa la aplicación
 */
function initApp() {
  console.log('Numerical Calculator JS - Starting...');
  
  // Inicializar Desmos si está disponible
  if (typeof Desmos !== 'undefined') {
    initDesmos();
  }
  
  // Configurar event listeners
  setupEventListeners();
  
  // Renderizar parámetros iniciales
  const checkedRadio = document.querySelector('input[name="method"]:checked');
  if (checkedRadio) {
    renderParameters(checkedRadio.value);
    selectMethodCard(checkedRadio.value);
  }
  
  console.log('Application initialized successfully');
}

/**
 * Configura los event listeners
 */
function setupEventListeners() {
  // Formulario de cálculo
  const form = document.querySelector('form');
  if (form) {
    form.addEventListener('submit', handleCalculate);
  }
  
  // Cambios en métodos
  document.querySelectorAll('input[name="method"]').forEach(radio => {
    radio.addEventListener('change', (e) => {
      AppState.currentMethod = e.target.value;
      renderParameters(e.target.value);
      updateGraph();
    });
  });
}

/**
 * Maneja el evento de cálculo
 */
async function handleCalculate(e) {
  e.preventDefault();
  
  try {
    const method = AppState.currentMethod;
    const formData = new FormData(e.target);
    
    // Mostrar indicador de carga
    showLoading();
    
    // Ejecutar método según tipo
    let result;
    
    if (isRootFindingMethod(method)) {
      result = await executeRootMethod(method, formData);
    } else if (isLinearSystemMethod(method)) {
      result = await executeLinearMethod(method, formData);
    } else if (isInterpolationMethod(method)) {
      result = await executeInterpolationMethod(method, formData);
    } else {
      throw new Error('Unknown method: ' + method);
    }
    
    // Guardar resultado
    AppState.lastResult = result;
    
    // Mostrar resultado
    displayResult(result, method);
    
    // Ocultar indicador de carga
    hideLoading();
    
  } catch (error) {
    hideLoading();
    displayError(error.message);
    console.error('Calculation error:', error);
  }
}

/**
 * Ejecuta un método de búsqueda de raíces
 */
async function executeRootMethod(method, formData) {
  const expr = formData.get('function');
  const f = createFunction(expr);
  
  let result;
  
  switch (method) {
    case 'incremental_search':
      const x0 = parseFloat(formData.get('x0'));
      const delta = parseFloat(formData.get('delta'));
      const nmax = parseInt(formData.get('nmax'));
      result = incrementalSearch(f, x0, delta, nmax);
      break;
      
    case 'bisection':
      const a = parseFloat(formData.get('a'));
      const b = parseFloat(formData.get('b'));
      const tol = parseFloat(formData.get('tolerance'));
      const nmax_b = parseInt(formData.get('nmax'));
      result = bisectionMethod(f, a, b, tol, nmax_b);
      break;
      
    case 'false_position':
      const a_fp = parseFloat(formData.get('a'));
      const b_fp = parseFloat(formData.get('b'));
      const tol_fp = parseFloat(formData.get('tolerance'));
      const nmax_fp = parseInt(formData.get('nmax'));
      result = falsePositionMethod(f, a_fp, b_fp, tol_fp, nmax_fp);
      break;
      
    case 'fixed_point':
      const gExpr = formData.get('g_function');
      const g = createFunction(gExpr);
      const x0_fp = parseFloat(formData.get('x0'));
      const tol_fpm = parseFloat(formData.get('tolerance'));
      const nmax_fpm = parseInt(formData.get('nmax'));
      result = fixedPointMethod(g, x0_fp, tol_fpm, nmax_fpm);
      result.f = f; // Para graficar
      break;
      
    case 'newton':
      const dfExpr = formData.get('derivative');
      const df = dfExpr ? createFunction(dfExpr) : null;
      const x0_n = parseFloat(formData.get('x0'));
      const tol_n = parseFloat(formData.get('tolerance'));
      const nmax_n = parseInt(formData.get('nmax'));
      result = newtonMethod(f, df, x0_n, tol_n, nmax_n);
      break;
      
    case 'secant':
      const x0_s = parseFloat(formData.get('x0'));
      const x1_s = parseFloat(formData.get('x1'));
      const tol_s = parseFloat(formData.get('tolerance'));
      const nmax_s = parseInt(formData.get('nmax'));
      result = secantMethod(f, x0_s, x1_s, tol_s, nmax_s);
      break;
      
    case 'newton_multiple':
      const dfExpr_m = formData.get('derivative');
      const d2fExpr_m = formData.get('second_derivative');
      const df_m = dfExpr_m ? createFunction(dfExpr_m) : null;
      const d2f_m = d2fExpr_m ? createFunction(d2fExpr_m) : null;
      const x0_nm = parseFloat(formData.get('x0'));
      const tol_nm = parseFloat(formData.get('tolerance'));
      const nmax_nm = parseInt(formData.get('nmax'));
      result = newtonMultipleRoots(f, df_m, d2f_m, x0_nm, tol_nm, nmax_nm);
      break;
  }
  
  result.method = method;
  result.f = f;
  result.expr = expr;
  
  return result;
}

/**
 * Ejecuta un método de sistemas lineales
 */
async function executeLinearMethod(method, formData) {
  let result;
  
  if (method === 'forward_substitution' || method === 'backward_substitution') {
    // Estos usan matriz aumentada
    const M = parseMatrix(formData.get('augmented_matrix'));
    
    if (method === 'forward_substitution') {
      result = forwardSubstitutionMethod(M);
    } else {
      result = backwardSubstitutionMethod(M);
    }
  } else {
    // Todos los demás usan A y b separados
    const A = parseMatrix(formData.get('matrix_A'));
    const b = parseVector(formData.get('vector_b'));
    
    switch (method) {
      case 'gauss':
        result = gaussElimination(A, b, true);
        break;
        
      case 'partial_pivot':
        result = partialPivoting(A, b, true);
        break;
        
      case 'total_pivot':
        result = totalPivoting(A, b, true);
        break;
        
      case 'lu_simple':
        result = luSimple(A, b);
        break;
        
      case 'lu_partial':
        result = luPartialPivot(A, b);
        break;
        
      case 'crout':
        result = crout(A, b);
        break;
        
      case 'doolittle':
        result = doolittle(A, b);
        break;
        
      case 'cholesky':
        result = cholesky(A, b);
        break;
        
      case 'jacobi':
        const x0_j = parseVector(formData.get('x0_vector') || '0 '.repeat(b.length));
        const tol_j = parseFloat(formData.get('tolerance'));
        const nmax_j = parseInt(formData.get('nmax'));
        result = jacobi(A, b, x0_j, tol_j, nmax_j);
        break;
        
      case 'gauss_seidel':
        const x0_gs = parseVector(formData.get('x0_vector') || '0 '.repeat(b.length));
        const tol_gs = parseFloat(formData.get('tolerance'));
        const nmax_gs = parseInt(formData.get('nmax'));
        result = gaussSeidel(A, b, x0_gs, tol_gs, nmax_gs);
        break;
        
      case 'sor':
        const x0_sor = parseVector(formData.get('x0_vector') || '0 '.repeat(b.length));
        const omega = parseFloat(formData.get('omega'));
        const tol_sor = parseFloat(formData.get('tolerance'));
        const nmax_sor = parseInt(formData.get('nmax'));
        result = sor(A, b, x0_sor, omega, tol_sor, nmax_sor);
        break;
    }
  }
  
  result.method = method;
  return result;
}

/**
 * Ejecuta un método de interpolación
 */
async function executeInterpolationMethod(method, formData) {
  const X = parseVector(formData.get('x_values'));
  const Y = parseVector(formData.get('y_values'));
  
  if (X.length !== Y.length) {
    throw new Error('X and Y must have the same number of points');
  }
  
  let result;
  
  switch (method) {
    case 'vandermonde':
      result = vandermonde(X, Y);
      break;
      
    case 'divided_differences':
      result = dividedDifferences(X, Y);
      break;
      
    case 'lagrange':
      result = lagrange(X, Y);
      break;
      
    case 'linear_spline':
      result = linearSpline(X, Y);
      break;
      
    case 'quadratic_spline':
      result = quadraticSpline(X, Y);
      break;
      
    case 'cubic_spline':
      result = cubicSpline(X, Y);
      break;
  }
  
  result.method = method;
  result.X = X;
  result.Y = Y;
  
  return result;
}

/**
 * Muestra el resultado en la interfaz
 */
function displayResult(result, method) {
  const outputDiv = document.querySelector('.output');
  if (!outputDiv) return;
  
  let html = '';
  
  if (isRootFindingMethod(method)) {
    html = formatRootResult(result, method);
  } else if (isLinearSystemMethod(method)) {
    html = formatLinearResult(result, method);
  } else if (isInterpolationMethod(method)) {
    html = formatInterpolationResult(result, method);
  }
  
  outputDiv.innerHTML = html;
}

/**
 * Muestra un error
 */
function displayError(message) {
  const outputDiv = document.querySelector('.output');
  if (outputDiv) {
    outputDiv.innerHTML = `<div style="color: #ef4444;">Error: ${message}</div>`;
  }
}

/**
 * Muestra indicador de carga
 */
function showLoading() {
  const outputDiv = document.querySelector('.output');
  if (outputDiv) {
    outputDiv.innerHTML = '<div style="color: #2196F3;">Calculating... ⏳</div>';
  }
}

/**
 * Oculta indicador de carga
 */
function hideLoading() {
  // Se llama después de displayResult, no hace nada
}

/**
 * Selecciona visualmente una tarjeta de método
 */
function selectMethodCard(method) {
  document.querySelectorAll('.method-card').forEach(card => {
    const radio = card.querySelector('input[type="radio"]');
    if (radio && radio.value === method) {
      card.classList.add('selected');
    } else {
      card.classList.remove('selected');
    }
  });
}

// Funciones de verificación de tipo de método
function isRootFindingMethod(method) {
  return ['incremental_search', 'bisection', 'false_position', 'fixed_point', 
          'newton', 'secant', 'newton_multiple'].includes(method);
}

function isLinearSystemMethod(method) {
  return ['gauss', 'partial_pivot', 'total_pivot', 'lu_simple', 'lu_partial',
          'crout', 'doolittle', 'cholesky', 'jacobi', 'gauss_seidel', 'sor',
          'forward_substitution', 'backward_substitution'].includes(method);
}

function isInterpolationMethod(method) {
  return ['vandermonde', 'divided_differences', 'lagrange', 'linear_spline',
          'quadratic_spline', 'cubic_spline'].includes(method);
}

// Inicializar cuando el DOM esté listo
if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', initApp);
} else {
  initApp();
}
