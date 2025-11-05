/**
 * INTERFAZ DE USUARIO
 * Maneja renderizado de parámetros, tabs, graficación con Desmos
 */

// Definición de parámetros para cada método
const methodParams = {
  // Ecuaciones - Métodos de Intervalo
  'incremental_search': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2', required: true},
    {name: 'x0', label: 'Inicial x₀', type: 'number', value: '0', step: 'any', required: true},
    {name: 'delta', label: 'Delta (Δ)', type: 'number', value: '0.5', step: 'any', required: true},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50', required: true}
  ],
  'bisection': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2', required: true},
    {name: 'a', label: 'Inicio intervalo (a)', type: 'number', value: '1', step: 'any', required: true},
    {name: 'b', label: 'Fin intervalo (b)', type: 'number', value: '2', step: 'any', required: true},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any', required: true},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50', required: true}
  ],
  'false_position': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2', required: true},
    {name: 'a', label: 'Inicio intervalo (a)', type: 'number', value: '1', step: 'any', required: true},
    {name: 'b', label: 'Fin intervalo (b)', type: 'number', value: '2', step: 'any', required: true},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any', required: true},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50', required: true}
  ],

  // Métodos Abiertos
  'fixed_point': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2'},
    {name: 'g_function', label: 'Función g(x)', type: 'text', placeholder: 'ej: (x^3 - 2)/x'},
    {name: 'x0', label: 'Inicial x₀', type: 'number', value: '1', step: 'any'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'}
  ],
  'newton': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2'},
    {name: 'x0', label: 'Inicial x₀', type: 'number', value: '1.5', step: 'any'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'},
    {name: 'derivative', label: 'Derivada f\'(x) (opcional)', type: 'text', placeholder: 'ej: 3*x^2 - 1', required: false}
  ],
  'secant': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2'},
    {name: 'x0', label: 'Primer punto x₀', type: 'number', value: '1', step: 'any'},
    {name: 'x1', label: 'Segundo punto x₁', type: 'number', value: '2', step: 'any'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'}
  ],
  'newton_multiple': [
    {name: 'function', label: 'Función f(x)', type: 'text', placeholder: 'ej: x^3 - x - 2', value: 'x^3 - x - 2'},
    {name: 'x0', label: 'Inicial x₀', type: 'number', value: '1.5', step: 'any'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'},
    {name: 'derivative', label: 'Primera derivada f\'(x) (opcional)', type: 'text', placeholder: 'ej: 3*x^2 - 1', required: false},
    {name: 'second_derivative', label: 'Segunda derivada f\'\'(x) (opcional)', type: 'text', placeholder: 'ej: 6*x', required: false}
  ],

  // Sistemas Lineales - Métodos Directos
  'gauss': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'partial_pivot': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'total_pivot': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'lu_simple': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'lu_partial': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'crout': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'doolittle': [
    {name: 'matrix_A', label: 'Matriz A (cada fila en una línea)', type: 'textarea', placeholder: '2 3 -1\n1 -2 3\n0 1 1', value: '2 3 -1\n1 -2 3\n0 1 1'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],
  'cholesky': [
    {name: 'matrix_A', label: 'Matriz A (simétrica def. positiva)', type: 'textarea', placeholder: '4 2 1\n2 5 3\n1 3 6', value: '4 2 1\n2 5 3\n1 3 6'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '1 -4 2'}
  ],

  // Sistemas Lineales - Métodos Iterativos
  'jacobi': [
    {name: 'matrix_A', label: 'Matriz A (diagonalmente dominante)', type: 'textarea', placeholder: '4 -1 0\n-1 4 -1\n0 -1 4', value: '4 -1 0\n-1 4 -1\n0 -1 4'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '15 10 10'},
    {name: 'x0_vector', label: 'Aproximación inicial x₀ (separado por espacios)', type: 'text', value: '0 0 0'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'}
  ],
  'gauss_seidel': [
    {name: 'matrix_A', label: 'Matriz A (diagonalmente dominante)', type: 'textarea', placeholder: '4 -1 0\n-1 4 -1\n0 -1 4', value: '4 -1 0\n-1 4 -1\n0 -1 4'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '15 10 10'},
    {name: 'x0_vector', label: 'Aproximación inicial x₀ (separado por espacios)', type: 'text', value: '0 0 0'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'}
  ],
  'sor': [
    {name: 'matrix_A', label: 'Matriz A (diagonalmente dominante)', type: 'textarea', placeholder: '4 -1 0\n-1 4 -1\n0 -1 4', value: '4 -1 0\n-1 4 -1\n0 -1 4'},
    {name: 'vector_b', label: 'Vector b (separado por espacios)', type: 'text', value: '15 10 10'},
    {name: 'x0_vector', label: 'Aproximación inicial x₀ (separado por espacios)', type: 'text', value: '0 0 0'},
    {name: 'omega', label: 'Factor de relajación ω (0 < ω < 2)', type: 'number', value: '1.5', step: '0.1', min: '0.01', max: '1.99'},
    {name: 'tolerance', label: 'Tolerancia', type: 'number', value: '0.000001', step: 'any'},
    {name: 'nmax', label: 'Iteraciones máx (Nmax)', type: 'number', value: '50'}
  ],

  // Interpolación
  'vandermonde': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],
  'divided_differences': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],
  'lagrange': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],
  'linear_spline': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],
  'quadratic_spline': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],
  'cubic_spline': [
    {name: 'x_values', label: 'Valores X (separados por espacios)', type: 'text', value: '0 1 2 3', placeholder: 'ej: 0 1 2 3'},
    {name: 'y_values', label: 'Valores Y (separados por espacios)', type: 'text', value: '1 2 0 4', placeholder: 'ej: 1 2 0 4'}
  ],

  // Métodos Adicionales
  'forward_substitution': [
    {name: 'augmented_matrix', label: 'Matriz aumentada [L|b] (triangular inferior, cada fila en una línea)', type: 'textarea', placeholder: '2 0 0 4\n3 5 0 19\n1 2 4 24', value: '2 0 0 4\n3 5 0 19\n1 2 4 24'}
  ],
  'backward_substitution': [
    {name: 'augmented_matrix', label: 'Matriz aumentada [U|b] (triangular superior, cada fila en una línea)', type: 'textarea', placeholder: '4 2 1 11\n0 5 3 19\n0 0 2 6', value: '4 2 1 11\n0 5 3 19\n0 0 2 6'}
  ]
};

/**
 * Renderiza los parámetros dinámicamente según el método
 */
function renderParameters(method) {
  const container = document.getElementById('dynamic-params');
  if (!container) return;
  
  container.innerHTML = '';
  
  const params = methodParams[method];
  if (!params) {
    container.innerHTML = '<p style="color:var(--text-light);">No se necesitan parámetros para este método.</p>';
    return;
  }
  
  params.forEach(param => {
    const div = document.createElement('div');
    div.style.marginBottom = '10px';
    
    const label = document.createElement('label');
    label.textContent = param.label;
    label.style.display = 'block';
    label.style.marginBottom = '4px';
    div.appendChild(label);
    
    let input;
    if (param.type === 'textarea') {
      input = document.createElement('textarea');
      input.rows = 3;
    } else {
      input = document.createElement('input');
      input.type = param.type;
    }
    
    input.name = param.name;
    if (param.placeholder) input.placeholder = param.placeholder;
    if (param.value !== undefined) input.value = param.value;
    if (param.step) input.step = param.step;
    if (param.min !== undefined) input.min = param.min;
    if (param.max !== undefined) input.max = param.max;
    if (param.required !== false) input.required = true;
    
    // Agregar listeners para actualizar gráfica
    if (param.name === 'function' || param.name === 'g_function' || 
        param.name === 'a' || param.name === 'b' || param.name === 'x0') {
      input.addEventListener('input', updateGraph);
      input.addEventListener('change', updateGraph);
    }
    
    div.appendChild(input);
    container.appendChild(div);
  });
  
  // Actualizar gráfica después de renderizar
  setTimeout(updateGraph, 50);
}

/**
 * Cambia entre tabs
 */
function switchTab(event, tabId) {
  document.querySelectorAll('.tab-button').forEach(btn => btn.classList.remove('active'));
  document.querySelectorAll('.tab-content').forEach(content => content.classList.remove('active'));
  event.currentTarget.classList.add('active');
  document.getElementById(tabId).classList.add('active');
}

/**
 * Selecciona un método
 */
function selectMethod(element) {
  document.querySelectorAll('.method-card').forEach(card => card.classList.remove('selected'));
  element.classList.add('selected');
  
  const radio = element.querySelector('input[type="radio"]');
  radio.checked = true;
  
  AppState.currentMethod = radio.value;
  renderParameters(radio.value);
}

// ============================================
// INTEGRACIÓN CON DESMOS
// ============================================

/**
 * Inicializa Desmos
 */
function initDesmos() {
  const calculatorElement = document.getElementById('calculator');
  if (calculatorElement && !AppState.calculator) {
    try {
      AppState.calculator = Desmos.GraphingCalculator(calculatorElement, {
        keypad: false,
        expressions: false,
        settingsMenu: false,
        zoomButtons: true,
        expressionsTopbar: false,
        border: false
      });
      console.log('Desmos calculator initialized');
    } catch (e) {
      console.warn('Could not initialize Desmos:', e);
    }
  }
}

/**
 * Convierte sintaxis Python/SymPy a LaTeX para Desmos
 */
function pythonToLatex(expr) {
  if (!expr) return '';
  
  let latex = expr.toString();
  
  // Reemplazar potencias
  latex = latex.replace(/\*\*/g, '^');
  latex = latex.replace(/\^/g, '^');
  
  // Reemplazar ln y log
  latex = latex.replace(/ln\(/g, '\\ln(');
  latex = latex.replace(/log\(/g, '\\log(');
  
  // Reemplazar funciones trigonométricas
  latex = latex.replace(/sin\(/g, '\\sin(');
  latex = latex.replace(/cos\(/g, '\\cos(');
  latex = latex.replace(/tan\(/g, '\\tan(');
  latex = latex.replace(/asin\(/g, '\\arcsin(');
  latex = latex.replace(/acos\(/g, '\\arccos(');
  latex = latex.replace(/atan\(/g, '\\arctan(');
  
  // Reemplazar funciones hiperbólicas
  latex = latex.replace(/sinh\(/g, '\\sinh(');
  latex = latex.replace(/cosh\(/g, '\\cosh(');
  latex = latex.replace(/tanh\(/g, '\\tanh(');
  
  // Reemplazar exp y sqrt
  latex = latex.replace(/exp\(/g, 'e^{(');
  latex = latex.replace(/sqrt\(/g, '\\sqrt{');
  
  // Limpiar multiplicaciones
  latex = latex.replace(/\*/g, ' \\cdot ');
  
  return latex;
}

/**
 * Actualiza la gráfica según la función de entrada
 */
function updateGraph() {
  if (!AppState.calculator) return;
  
  const functionInput = document.querySelector('input[name="function"]');
  const gFunctionInput = document.querySelector('input[name="g_function"]');
  
  if (!functionInput) {
    document.getElementById('graph-section').style.display = 'none';
    return;
  }
  
  const funcValue = functionInput.value.trim();
  
  if (!funcValue) {
    document.getElementById('graph-section').style.display = 'none';
    return;
  }
  
  // Mostrar sección de gráfica
  document.getElementById('graph-section').style.display = 'block';
  
  // Limpiar expresiones previas
  AppState.calculator.setBlank();
  
  // Convertir y graficar f(x)
  const latexFunc = pythonToLatex(funcValue);
  AppState.calculator.setExpression({
    id: 'main-function',
    latex: `y=${latexFunc}`,
    color: '#2563eb',
    lineWidth: 2
  });
  
  // Agregar eje x
  AppState.calculator.setExpression({
    id: 'x-axis',
    latex: 'y=0',
    color: '#9ca3af',
    lineStyle: Desmos.Styles.DASHED,
    lineWidth: 1
  });
  
  // Graficar g(x) si existe (para punto fijo)
  if (gFunctionInput && gFunctionInput.value.trim()) {
    const gLatex = pythonToLatex(gFunctionInput.value.trim());
    AppState.calculator.setExpression({
      id: 'g-function',
      latex: `y=${gLatex}`,
      color: '#16a34a',
      lineWidth: 2
    });
    
    // Agregar línea y=x para punto fijo
    AppState.calculator.setExpression({
      id: 'identity',
      latex: 'y=x',
      color: '#dc2626',
      lineStyle: Desmos.Styles.DASHED,
      lineWidth: 1
    });
  }
  
  // Agregar marcadores de intervalo para bisección/regla falsa
  const aInput = document.querySelector('input[name="a"]');
  const bInput = document.querySelector('input[name="b"]');
  
  if (aInput && bInput && aInput.value && bInput.value) {
    const a = parseFloat(aInput.value);
    const b = parseFloat(bInput.value);
    
    if (!isNaN(a) && !isNaN(b)) {
      AppState.calculator.setExpression({
        id: 'interval-a',
        latex: `x=${a}`,
        color: '#f59e0b',
        lineStyle: Desmos.Styles.DASHED,
        lineWidth: 1
      });
      
      AppState.calculator.setExpression({
        id: 'interval-b',
        latex: `x=${b}`,
        color: '#f59e0b',
        lineStyle: Desmos.Styles.DASHED,
        lineWidth: 1
      });
    }
  }
  
  // Agregar punto inicial para Newton/Secante
  const x0Input = document.querySelector('input[name="x0"]');
  if (x0Input && x0Input.value) {
    const x0 = parseFloat(x0Input.value);
    if (!isNaN(x0)) {
      AppState.calculator.setExpression({
        id: 'initial-point',
        latex: `x=${x0}`,
        color: '#8b5cf6',
        lineStyle: Desmos.Styles.DASHED,
        lineWidth: 1
      });
    }
  }
}
