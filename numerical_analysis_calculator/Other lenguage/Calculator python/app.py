from flask import Flask, render_template, request, redirect, url_for
import numpy as np
import sympy as sp
from sympy import symbols, lambdify, diff
import math

app = Flask(__name__)

def create_function(expr):
    """Crea una función evaluable a partir de una expresión string"""
    x = symbols('x')
    expr = expr.replace('ln(', 'log(')
    expr = expr.replace('Log(', 'log(')
    parsed = sp.sympify(expr)
    return lambdify(x, parsed, 'numpy')

def numeric_derivative(f, x, h=1e-6):
    """Derivada numérica usando diferencias centradas"""
    return (f(x + h) - f(x - h)) / (2 * h)

def numeric_second_derivative(f, x, h=1e-5):
    """Segunda derivada numérica"""
    return (f(x + h) - 2*f(x) + f(x - h)) / (h**2)

# MÉTODOS NUMÉRICOS

def incremental_search(f, x0, delta, nmax):
    """Búsqueda incremental"""
    iterations = []
    fx0 = f(x0)

    for i in range(nmax):
        x1 = x0 + delta
        fx1 = f(x1)
        iterations.append([i, x0, fx0, x1, fx1])

        if fx0 * fx1 < 0:
            return {"interval": [x0, x1], "iterations": iterations}

        x0, fx0 = x1, fx1

    return {"interval": None, "iterations": iterations}

def bisection_method(f, a, b, tol, nmax):
    """Método de bisección"""
    iterations = []

    for i in range(nmax):
        c = (a + b) / 2
        fc = f(c)
        iterations.append([i, a, b, c, fc])

        if abs(fc) < tol or (b - a) / 2 < tol:
            return {"root": c, "iterations": iterations}

        if f(a) * fc < 0:
            b = c
        else:
            a = c

    return {"root": (a + b) / 2, "iterations": iterations}

def false_position_method(f, a, b, tol, nmax):
    """Método de falsa posición"""
    iterations = []

    for i in range(nmax):
        fa, fb = f(a), f(b)
        c = (a * fb - b * fa) / (fb - fa)
        fc = f(c)
        iterations.append([i, a, b, c, fc])

        if abs(fc) < tol:
            return {"root": c, "iterations": iterations}

        if fa * fc < 0:
            b = c
        else:
            a = c

    return {"root": c, "iterations": iterations}

def fixed_point_method(g, x0, tol, nmax):
    """Método de punto fijo"""
    iterations = []

    for i in range(nmax):
        x1 = g(x0)
        iterations.append([i, x0, x1])

        if abs(x1 - x0) < tol:
            return {"root": x1, "iterations": iterations}

        x0 = x1

    return {"root": x0, "iterations": iterations}

def newton_method(f, df, x0, tol, nmax):
    """Método de Newton-Raphson"""
    iterations = []

    for i in range(nmax):
        fx = f(x0)
        dfx = df(x0) if df else numeric_derivative(f, x0)

        if dfx == 0:
            break

        x1 = x0 - fx / dfx
        iterations.append([i, x0, fx, dfx, x1])

        if abs(x1 - x0) < tol:
            return {"root": x1, "iterations": iterations}

        x0 = x1

    return {"root": x0, "iterations": iterations}

def secant_method(f, x0, x1, tol, nmax):
    """Método de la secante"""
    iterations = []

    for i in range(nmax):
        fx0, fx1 = f(x0), f(x1)

        if fx1 - fx0 == 0:
            break

        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        iterations.append([i, x0, x1, x2, fx1])

        if abs(x2 - x1) < tol:
            return {"root": x2, "iterations": iterations}

        x0, x1 = x1, x2

    return {"root": x1, "iterations": iterations}

def newton_multiple_roots(f, df, d2f, x0, tol, nmax):
    """Método de Newton para raíces múltiples"""
    iterations = []

    for i in range(nmax):
        fx = f(x0)
        dfx = df(x0) if df else numeric_derivative(f, x0)
        d2fx = d2f(x0) if d2f else numeric_second_derivative(f, x0)

        denom = dfx**2 - fx * d2fx
        if denom == 0:
            break

        x1 = x0 - (fx * dfx) / denom
        iterations.append([i, x0, fx, dfx, d2fx, x1])

        if abs(x1 - x0) < tol:
            return {"root": x1, "iterations": iterations}

        x0 = x1

    return {"root": x0, "iterations": iterations}

def parse_matrix(matrix_text):
    """Parsea una matriz desde texto"""
    rows = [row.strip() for row in matrix_text.strip().split('\n') if row.strip()]
    return [[float(x) for x in row.split()] for row in rows]

def gauss_elimination(A, b):
    """Eliminación Gaussiana"""
    n = len(A)
    # Crear matriz aumentada
    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    # Eliminación hacia adelante
    for k in range(n - 1):
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

    # Sustitución hacia atrás
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    return x

def partial_pivoting(A, b):
    """Eliminación Gaussiana con pivoteo parcial"""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    for k in range(n - 1):
        # Encontrar fila con el mayor elemento en columna k
        max_row = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > abs(M[max_row][k]):
                max_row = i

        # Intercambiar filas
        if max_row != k:
            M[k], M[max_row] = M[max_row], M[k]

        # Eliminación
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

    # Sustitución hacia atrás
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    return x

def total_pivoting(A, b):
    """Eliminación Gaussiana con pivoteo total"""
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    marks = list(range(n))

    for k in range(n - 1):
        # Encontrar el elemento máximo en la submatriz
        max_val, max_row, max_col = 0, k, k
        for i in range(k, n):
            for j in range(k, n):
                if abs(M[i][j]) > max_val:
                    max_val, max_row, max_col = abs(M[i][j]), i, j

        if max_val == 0:
            continue

        # Intercambiar filas
        if max_row != k:
            M[k], M[max_row] = M[max_row], M[k]

        # Intercambiar columnas
        if max_col != k:
            for i in range(n):
                M[i][k], M[i][max_col] = M[i][max_col], M[i][k]
            marks[k], marks[max_col] = marks[max_col], marks[k]

        # Eliminación
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

    # Sustitución hacia atrás
    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    # Reordenar solución
    x_final = [0] * n
    for i in range(n):
        x_final[marks[i]] = x[i]

    return x_final

@app.route('/')
def index():
    return render_template('index.html',
                         result=None,
                         iterations=None,
                         headers=None,
                         error=None,
                         method='bisection')

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        method = request.form['method']
        result = None
        iterations = []
        headers = []

        if method in ['bisection', 'false_position', 'incremental_search', 'fixed_point', 'newton', 'secant', 'newton_multiple']:
            # Métodos para encontrar raíces
            expr = request.form['function']
            f = create_function(expr)

            if method == 'incremental_search':
                x0 = float(request.form['x0'])
                delta = float(request.form['delta'])
                nmax = int(request.form['nmax'])
                res = incremental_search(f, x0, delta, nmax)

                result = f"Intervalo con cambio de signo: {res['interval']}" if res['interval'] else "No se encontró intervalo en Nmax iteraciones"
                iterations = res['iterations']
                headers = ['i', 'x0', 'f(x0)', 'x1', 'f(x1)']

            elif method == 'bisection':
                a, b = float(request.form['a']), float(request.form['b'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = bisection_method(f, a, b, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'a', 'b', 'c', 'f(c)']

            elif method == 'false_position':
                a, b = float(request.form['a']), float(request.form['b'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = false_position_method(f, a, b, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'a', 'b', 'c', 'f(c)']

            elif method == 'fixed_point':
                g_expr = request.form['g_function']
                g = create_function(g_expr)
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = fixed_point_method(g, x0, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'x_old', 'x_new']

            elif method == 'newton':
                df = create_function(request.form['derivative']) if request.form.get('derivative') else None
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = newton_method(f, df, x0, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'x', 'f(x)', "f'(x)", 'x_next']

            elif method == 'secant':
                x0, x1 = float(request.form['x0']), float(request.form['x1'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = secant_method(f, x0, x1, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'x0', 'x1', 'x2', 'f(x1)']

            elif method == 'newton_multiple':
                df = create_function(request.form['derivative']) if request.form.get('derivative') else None
                d2f = create_function(request.form['second_derivative']) if request.form.get('second_derivative') else None
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = newton_multiple_roots(f, df, d2f, x0, tol, nmax)

                result = f"Raíz aproximada: {res['root']}"
                iterations = res['iterations']
                headers = ['i', 'x', 'f(x)', "f'", "f''", 'x_next']

        elif method in ['gauss', 'partial_pivot', 'total_pivot']:
            # Métodos para sistemas lineales
            A = parse_matrix(request.form['matrix_A'])
            b = [float(x) for x in request.form['vector_b'].split()]

            if method == 'gauss':
                solution = gauss_elimination(A, b)
            elif method == 'partial_pivot':
                solution = partial_pivoting(A, b)
            elif method == 'total_pivot':
                solution = total_pivoting(A, b)

            result = f"Solución: {', '.join([f'{x:.8g}' for x in solution])}"

        return render_template('index.html',
                             result=result,
                             iterations=iterations,
                             headers=headers,
                             error=None,
                             method=method)

    except Exception as e:
        return render_template('index.html',
                             result=None,
                             iterations=None,
                             headers=None,
                             error=str(e),
                             method=request.form.get('method', 'bisection'))

if __name__ == '__main__':
    app.run(debug=True)