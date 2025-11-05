from flask import Flask, render_template, request, redirect, url_for
import numpy as np
import sympy as sp
from sympy import symbols, lambdify, diff
import math

app = Flask(__name__)

def create_function(expr):
    x = symbols('x')
    expr = expr.replace('ln(', 'log(')
    expr = expr.replace('Log(', 'log(')
    parsed = sp.sympify(expr)
    return lambdify(x, parsed, 'numpy')

def numeric_derivative(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2 * h)

def numeric_second_derivative(f, x, h=1e-5):
    return (f(x + h) - 2*f(x) + f(x - h)) / (h**2)

def incremental_search(f, x0, delta, nmax):
    iterations = []
    intervals = []

    xant = x0
    fant = f(xant)
    xact = xant + delta
    fact = f(xact)

    for i in range(nmax):
        iterations.append([i, xant, fant, xact, fact])

        if fant * fact < 0:
            intervals.append([xant, xact])

        xant = xact
        fant = fact
        xact = xant + delta
        fact = f(xact)

    return {"intervals": intervals, "iterations": iterations}

def bisection_method(f, a, b, tol, nmax):
    iterations = []

    fa = f(a)
    pm = (a + b) / 2
    fpm = f(pm)
    E = 1000
    cont = 1

    iterations.append([0, a, b, pm, fpm])

    while E > tol and cont < nmax:
        if fa * fpm < 0:
            b = pm
        else:
            a = pm

        p0 = pm
        pm = (a + b) / 2
        fpm = f(pm)
        E = abs(pm - p0)

        iterations.append([cont, a, b, pm, fpm])
        cont = cont + 1
        fa = f(a)

    return {"root": pm, "iterations": iterations, "error": E, "iter_count": cont}

def false_position_method(f, a, b, tol, nmax):
    iterations = []

    fa = f(a)
    fb = f(b)
    pm = (fb * a - fa * b) / (fb - fa)
    fpm = f(pm)
    E = 1000
    cont = 1

    iterations.append([0, a, b, pm, fpm])

    while E > tol and cont < nmax:
        if fa * fpm < 0:
            b = pm
        else:
            a = pm

        p0 = pm
        fa = f(a)
        fb = f(b)
        pm = (fb * a - fa * b) / (fb - fa)
        fpm = f(pm)
        E = abs(pm - p0)

        iterations.append([cont, a, b, pm, fpm])
        cont = cont + 1

    return {"root": pm, "iterations": iterations, "error": E, "iter_count": cont}

def fixed_point_method(g, x0, tol, nmax):
    iterations = []

    for i in range(nmax):
        x1 = g(x0)
        iterations.append([i, x0, x1])

        if abs(x1 - x0) < tol:
            return {"root": x1, "iterations": iterations}

        x0 = x1

    return {"root": x0, "iterations": iterations}

def newton_method(f, df, x0, tol, nmax):
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
    rows = [row.strip() for row in matrix_text.strip().split('\n') if row.strip()]
    return [[float(x) for x in row.split()] for row in rows]

def gauss_elimination(A, b, return_stages=False):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []

    if return_stages:
        stages.append(("Initial augmented matrix", [row[:] for row in M]))

    for k in range(n - 1):
        if abs(M[k][k]) < 1e-10:
            raise ValueError(f"Zero pivot encountered at position ({k},{k}). Matrix may be singular. Consider using pivoting methods.")

        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

        if return_stages:
            stages.append((f"After elimination step {k+1}", [row[:] for row in M]))

    x = [0] * n
    for i in range(n - 1, -1, -1):
        if abs(M[i][i]) < 1e-10:
            raise ValueError(f"Zero diagonal element at position ({i},{i}). Matrix may be singular.")

        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    if return_stages:
        return {"solution": x, "stages": stages}
    return x

def partial_pivoting(A, b, return_stages=False):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []

    if return_stages:
        stages.append(("Initial augmented matrix", [row[:] for row in M]))

    for k in range(n - 1):
        max_row = k
        for i in range(k + 1, n):
            if abs(M[i][k]) > abs(M[max_row][k]):
                max_row = i

        if max_row != k:
            M[k], M[max_row] = M[max_row], M[k]
            if return_stages:
                stages.append((f"After row swap (pivot): row {k} ↔ row {max_row}", [row[:] for row in M]))
        if abs(M[k][k]) < 1e-10:
            raise ValueError(f"Zero pivot encountered at position ({k},{k}). Matrix may be singular. Consider using total pivoting.")

        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

        if return_stages:
            stages.append((f"After elimination step {k+1}", [row[:] for row in M]))

    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    if return_stages:
        return {"solution": x, "stages": stages}
    return x

def total_pivoting(A, b, return_stages=False):
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    marks = list(range(n))

    stages = []

    if return_stages:
        stages.append(("Initial augmented matrix", [row[:] for row in M]))

    for k in range(n - 1):
        max_val, max_row, max_col = 0, k, k
        for i in range(k, n):
            for j in range(k, n):
                if abs(M[i][j]) > max_val:
                    max_val, max_row, max_col = abs(M[i][j]), i, j
        if max_val < 1e-10:
            raise ValueError(f"Zero pivot encountered at position ({k},{k}). Matrix is singular and cannot be solved.")

        if max_row != k:
            M[k], M[max_row] = M[max_row], M[k]
            if return_stages:
                stages.append((f"After row swap: row {k} ↔ row {max_row}", [row[:] for row in M]))

        if max_col != k:
            for i in range(n):
                M[i][k], M[i][max_col] = M[i][max_col], M[i][k]
            marks[k], marks[max_col] = marks[max_col], marks[k]
            if return_stages:
                stages.append((f"After column swap: col {k} ↔ col {max_col}", [row[:] for row in M]))

        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

        if return_stages:
            stages.append((f"After elimination step {k+1}", [row[:] for row in M]))

    x = [0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    x_final = [0] * n
    for i in range(n):
        x_final[marks[i]] = x[i]

    if return_stages:
        return {"solution": x_final, "stages": stages, "marks": marks}
    return x_final

def forward_substitution(L, b):
    n = len(L)
    x = [0] * n

    for i in range(n):
        s = b[i]
        for j in range(i):
            s -= L[i][j] * x[j]
        x[i] = s / L[i][i]

    return x

def backward_substitution(U, b):
    n = len(U)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        s = b[i]
        for j in range(i + 1, n):
            s -= U[i][j] * x[j]
        x[i] = s / U[i][i]

    return x

def forward_substitution_method(M):
    n = len(M)
    x = [0.0] * n
    stages = []

    stages.append(("Initial augmented matrix [L|b]", [row[:] for row in M]))

    x[0] = M[0][n] / M[0][0]

    iteration_details = []
    iteration_details.append(f"x[0] = {M[0][n]:.6f} / {M[0][0]:.6f} = {x[0]:.6f}")

    for i in range(1, n):
        s = M[i][n]
        for j in range(i):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

        iteration_details.append(f"x[{i}] = ({M[i][n]:.6f} - sum) / {M[i][i]:.6f} = {x[i]:.6f}")

    stages.append(("Forward substitution steps", iteration_details))

    return {"solution": x, "stages": stages}

def backward_substitution_method(M):
    n = len(M)
    x = [0.0] * n
    stages = []

    stages.append(("Initial augmented matrix [U|b]", [row[:] for row in M]))

    x[n-1] = M[n-1][n] / M[n-1][n-1]

    iteration_details = []
    iteration_details.append(f"x[{n-1}] = {M[n-1][n]:.6f} / {M[n-1][n-1]:.6f} = {x[n-1]:.6f}")

    for i in range(n-2, -1, -1):
        s = M[i][n]
        for j in range(i+1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

        iteration_details.append(f"x[{i}] = ({M[i][n]:.6f} - sum) / {M[i][i]:.6f} = {x[i]:.6f}")

    stages.append(("Backward substitution steps", iteration_details))

    return {"solution": x, "stages": stages}

def lu_simple(A, b):
    n = len(A)

    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0 for j in range(n)] for i in range(n)]

    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []
    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    for i in range(n - 1):
        for j in range(i + 1, n):
            if M[j][i] != 0:
                L[j][i] = M[j][i] / M[i][i]
                for k in range(i, n + 1):
                    M[j][k] -= L[j][i] * M[i][k]

        for k in range(i, n):
            U[i][k] = M[i][k]

        L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
        U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
        stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    for k in range(n):
        U[n-1][k] = M[n-1][k]

    z = forward_substitution(L, b)
    x = backward_substitution(U, z)

    return {"solution": x, "L": L, "U": U, "stages": stages}

def lu_partial_pivot(A, b):
    n = len(A)

    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[0.0 for j in range(n)] for i in range(n)]
    P = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []
    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    P_copy = [[P[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy, "P": P_copy})

    for i in range(n - 1):
        max_row = i
        max_val = abs(M[i][i])

        for k in range(i + 1, n):
            if abs(M[k][i]) > max_val:
                max_val = abs(M[k][i])
                max_row = k

        if max_row != i:
            M[i], M[max_row] = M[max_row], M[i]
            P[i], P[max_row] = P[max_row], P[i]

            if i > 0:
                for k in range(i):
                    L[i][k], L[max_row][k] = L[max_row][k], L[i][k]

        for j in range(i + 1, n):
            if M[j][i] != 0:
                L[j][i] = M[j][i] / M[i][i]
                for k in range(i, n + 1):
                    M[j][k] -= L[j][i] * M[i][k]

        for k in range(i, n):
            U[i][k] = M[i][k]

        L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
        U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
        P_copy = [[P[i][j] for j in range(n)] for i in range(n)]
        stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy, "P": P_copy})

    for k in range(n):
        U[n-1][k] = M[n-1][k]

    Pb = [sum(P[i][j] * b[j] for j in range(n)) for i in range(n)]

    z = forward_substitution(L, Pb)
    x = backward_substitution(U, z)

    return {"solution": x, "L": L, "U": U, "P": P, "stages": stages}

def crout(A, b):
    n = len(A)

    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []
    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    for i in range(n - 1):
        for j in range(i, n):
            sum_val = sum(L[j][k] * U[k][i] for k in range(i))
            L[j][i] = A[j][i] - sum_val

        for j in range(i + 1, n):
            sum_val = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = (A[i][j] - sum_val) / L[i][i]

        L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
        U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
        stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    sum_val = sum(L[n-1][k] * U[k][n-1] for k in range(n - 1))
    L[n-1][n-1] = A[n-1][n-1] - sum_val

    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    z = forward_substitution(L, b)
    x = backward_substitution(U, z)

    return {"solution": x, "L": L, "U": U, "stages": stages}

def doolittle(A, b):
    n = len(A)

    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []
    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    for i in range(n - 1):
        for j in range(i, n):
            sum_val = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = A[i][j] - sum_val

        for j in range(i + 1, n):
            sum_val = sum(L[j][k] * U[k][i] for k in range(i))
            L[j][i] = (A[j][i] - sum_val) / U[i][i]

        L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
        U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
        stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    sum_val = sum(L[n-1][k] * U[k][n-1] for k in range(n - 1))
    U[n-1][n-1] = A[n-1][n-1] - sum_val

    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    z = forward_substitution(L, b)
    x = backward_substitution(U, z)

    return {"solution": x, "L": L, "U": U, "stages": stages}

def cholesky(A, b):
    import math

    n = len(A)
    L = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    U = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

    M = [row[:] + [b[i]] for i, row in enumerate(A)]

    stages = []
    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    for i in range(n - 1):
        sum_val = sum(L[i][k] * U[k][i] for k in range(i))
        value = A[i][i] - sum_val

        if value < 0:
            raise ValueError(f"Matrix is not positive definite (negative value {value} at position ({i},{i}))")

        L[i][i] = math.sqrt(value)
        U[i][i] = L[i][i]

        for j in range(i + 1, n):
            sum_val = sum(L[j][k] * U[k][i] for k in range(i))
            L[j][i] = (A[j][i] - sum_val) / U[i][i]

        for j in range(i + 1, n):
            sum_val = sum(L[i][k] * U[k][j] for k in range(i))
            U[i][j] = (A[i][j] - sum_val) / L[i][i]

        L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
        U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
        stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    sum_val = sum(L[n-1][k] * U[k][n-1] for k in range(n - 1))
    value = A[n-1][n-1] - sum_val

    if value < 0:
        raise ValueError(f"Matrix is not positive definite (negative value {value} at position ({n-1},{n-1}))")

    L[n-1][n-1] = math.sqrt(value)
    U[n-1][n-1] = L[n-1][n-1]

    L_copy = [[L[i][j] for j in range(n)] for i in range(n)]
    U_copy = [[U[i][j] for j in range(n)] for i in range(n)]
    stages.append({"M": [row[:] for row in M], "L": L_copy, "U": U_copy})

    z = forward_substitution(L, b)
    x = backward_substitution(U, z)

    return {"solution": x, "L": L, "U": U, "stages": stages}

def jacobi(A, b, x0, tol, nmax):
    import numpy as np

    n = len(A)
    A_np = np.array(A, dtype=float)
    b_np = np.array(b, dtype=float)
    for i in range(n):
        if A[i][i] == 0:
            raise ValueError(f"Diagonal element A[{i}][{i}] is zero. Jacobi method cannot proceed.")

    D = np.diag(np.diag(A_np))
    L = -np.tril(A_np) + D
    U = -np.triu(A_np) + D
    T = np.linalg.inv(D) @ (L + U)
    C = np.linalg.inv(D) @ b_np

    eigenvalues = np.linalg.eigvals(T)
    spectral_radius = max(abs(ev) for ev in eigenvalues)

    xant = np.array(x0, dtype=float)
    E = 1000
    cont = 0

    iterations = []

    while E > tol and cont < nmax:
        xact = T @ xant + C
        E = np.linalg.norm(xant - xact)
        xant = xact
        cont = cont + 1

        iterations.append([cont-1, E] + xact.tolist())

    x = xact
    iter_count = cont
    err = E

    return {
        "solution": x.tolist(),
        "iterations": iterations,
        "converged": E < tol,
        "T": T.tolist(),
        "C": C.tolist(),
        "spectral_radius": spectral_radius
    }

def gauss_seidel(A, b, x0, tol, nmax):
    import numpy as np

    n = len(A)
    A_np = np.array(A, dtype=float)
    b_np = np.array(b, dtype=float)
    for i in range(n):
        if A[i][i] == 0:
            raise ValueError(f"Diagonal element A[{i}][{i}] is zero. Gauss-Seidel method cannot proceed.")

    D = np.diag(np.diag(A_np))
    L = -np.tril(A_np) + D
    U = -np.triu(A_np) + D
    T = np.linalg.inv(D - L) @ U
    C = np.linalg.inv(D - L) @ b_np

    eigenvalues = np.linalg.eigvals(T)
    spectral_radius = max(abs(ev) for ev in eigenvalues)

    xant = np.array(x0, dtype=float)
    E = 1000
    cont = 0

    iterations = []

    while E > tol and cont < nmax:
        xact = T @ xant + C
        E = np.linalg.norm(xant - xact)
        xant = xact
        cont = cont + 1

        iterations.append([cont-1, E] + xact.tolist())

    x = xact
    iter_count = cont
    err = E

    return {
        "solution": x.tolist(),
        "iterations": iterations,
        "converged": E < tol,
        "T": T.tolist(),
        "C": C.tolist(),
        "spectral_radius": spectral_radius
    }

def sor(A, b, x0, omega, tol, nmax):
    import numpy as np

    n = len(A)
    A_np = np.array(A, dtype=float)
    b_np = np.array(b, dtype=float)
    for i in range(n):
        if A[i][i] == 0:
            raise ValueError(f"Diagonal element A[{i}][{i}] is zero. SOR method cannot proceed.")
    if omega <= 0 or omega >= 2:
        raise ValueError(f"Relaxation factor omega={omega} is out of range. Must be 0 < omega < 2.")

    D = np.diag(np.diag(A_np))
    L = -np.tril(A_np) + D
    U = -np.triu(A_np) + D
    T = np.linalg.inv(D - omega * L) @ ((1 - omega) * D + omega * U)
    C = omega * np.linalg.inv(D - omega * L) @ b_np

    eigenvalues = np.linalg.eigvals(T)
    spectral_radius = max(abs(ev) for ev in eigenvalues)

    xant = np.array(x0, dtype=float)
    E = 1000
    cont = 0

    iterations = []

    while E > tol and cont < nmax:
        xact = T @ xant + C
        E = np.linalg.norm(xant - xact)
        xant = xact
        cont = cont + 1

        iterations.append([cont-1, E] + xact.tolist())

    x = xact
    iter_count = cont
    err = E

    return {
        "solution": x.tolist(),
        "iterations": iterations,
        "converged": E < tol,
        "T": T.tolist(),
        "C": C.tolist(),
        "spectral_radius": spectral_radius,
        "omega": omega
    }

def vandermonde(X, Y):
    n = len(X)

    A = []
    for i in range(n):
        row = [X[i]**(n - 1 - j) for j in range(n)]
        A.append(row)

    coef = total_pivoting(A, Y, return_stages=False)

    return {"coefficients": coef, "degree": n - 1}

def divided_differences(X, Y):
    n = len(X)

    D = [[0.0 for _ in range(n)] for _ in range(n)]

    for i in range(n):
        D[i][0] = Y[i]

    for j in range(1, n):
        for i in range(j, n):
            D[i][j] = (D[i][j-1] - D[i-1][j-1]) / (X[i] - X[i-j])

    coef = [D[i][i] for i in range(n)]

    return {"coefficients": coef, "x_points": X, "degree": n - 1, "table": D}

def lagrange(X, Y):
    n = len(X)
    L = np.zeros((n, n))

    for i in range(n):
        aux0 = [X[j] for j in range(n) if j != i]
        aux = np.array([1.0, -aux0[0]])

        for j in range(1, n - 1):
            aux = np.convolve(aux, np.array([1.0, -aux0[j]]))

        denominator = np.polyval(aux, X[i])
        L[i, :] = aux / denominator

    Coef = np.dot(Y, L)

    return {
        "coefficients": Coef,
        "lagrange_basis": L,
        "x_points": X,
        "y_points": Y,
        "degree": n - 1
    }

def linear_spline(X, Y):
    n = len(X)

    Coef = []
    for i in range(n - 1):
        if X[i + 1] == X[i]:
            raise ValueError(f"Duplicate x-values at index {i} and {i+1}: X[{i}] = X[{i+1}] = {X[i]}. X values must be distinct.")

        a_i = (Y[i + 1] - Y[i]) / (X[i + 1] - X[i])
        b_i = Y[i] - a_i * X[i]

        Coef.append([a_i, b_i])

    return {
        "coefficients": Coef,
        "x_points": X,
        "y_points": Y,
        "num_segments": n - 1
    }

def quadratic_spline(X, Y):
    n = len(X)
    m = 3 * (n - 1)
    A = [[0.0 for _ in range(m)] for _ in range(m)]
    b = [0.0 for _ in range(m)]

    for i in range(n - 1):
        row_idx = i + 1
        col_start = 3 * i
        A[row_idx][col_start] = X[i + 1]**2
        A[row_idx][col_start + 1] = X[i + 1]
        A[row_idx][col_start + 2] = 1
        b[row_idx] = Y[i + 1]

    A[0][0] = X[0]**2
    A[0][1] = X[0]
    A[0][2] = 1
    b[0] = Y[0]

    for i in range(1, n - 1):
        row_idx = n - 1 + i
        col_start = 3 * i - 3
        A[row_idx][col_start] = X[i]**2
        A[row_idx][col_start + 1] = X[i]
        A[row_idx][col_start + 2] = 1
        A[row_idx][col_start + 3] = -X[i]**2
        A[row_idx][col_start + 4] = -X[i]
        A[row_idx][col_start + 5] = -1
        b[row_idx] = 0

    for i in range(1, n - 1):
        row_idx = 2 * n - 3 + i
        col_start = 3 * i - 3
        A[row_idx][col_start] = 2 * X[i]
        A[row_idx][col_start + 1] = 1
        A[row_idx][col_start + 2] = 0
        A[row_idx][col_start + 3] = -2 * X[i]
        A[row_idx][col_start + 4] = -1
        A[row_idx][col_start + 5] = 0
        b[row_idx] = 0

    A[m - 1][0] = 2
    b[m - 1] = 0

    Saux = partial_pivoting(A, b, return_stages=False)

    Coef = []
    for i in range(n - 1):
        Coef.append([Saux[3 * i], Saux[3 * i + 1], Saux[3 * i + 2]])

    return {
        "coefficients": Coef,
        "x_points": X,
        "y_points": Y,
        "num_segments": n - 1
    }

def cubic_spline(X, Y):
    n = len(X)
    m = 4 * (n - 1)
    A = [[0.0 for _ in range(m)] for _ in range(m)]
    b = [0.0 for _ in range(m)]

    for i in range(n - 1):
        row_idx = i + 1
        col_start = 4 * i
        A[row_idx][col_start] = X[i + 1]**3
        A[row_idx][col_start + 1] = X[i + 1]**2
        A[row_idx][col_start + 2] = X[i + 1]
        A[row_idx][col_start + 3] = 1
        b[row_idx] = Y[i + 1]

    A[0][0] = X[0]**3
    A[0][1] = X[0]**2
    A[0][2] = X[0]
    A[0][3] = 1
    b[0] = Y[0]

    for i in range(1, n - 1):
        row_idx = n - 1 + i
        col_start = 4 * i - 4
        A[row_idx][col_start] = X[i]**3
        A[row_idx][col_start + 1] = X[i]**2
        A[row_idx][col_start + 2] = X[i]
        A[row_idx][col_start + 3] = 1
        A[row_idx][col_start + 4] = -X[i]**3
        A[row_idx][col_start + 5] = -X[i]**2
        A[row_idx][col_start + 6] = -X[i]
        A[row_idx][col_start + 7] = -1
        b[row_idx] = 0

    for i in range(1, n - 1):
        row_idx = 2 * n - 3 + i
        col_start = 4 * i - 4
        A[row_idx][col_start] = 3 * X[i]**2
        A[row_idx][col_start + 1] = 2 * X[i]
        A[row_idx][col_start + 2] = 1
        A[row_idx][col_start + 3] = 0
        A[row_idx][col_start + 4] = -3 * X[i]**2
        A[row_idx][col_start + 5] = -2 * X[i]
        A[row_idx][col_start + 6] = -1
        A[row_idx][col_start + 7] = 0
        b[row_idx] = 0

    for i in range(1, n - 1):
        row_idx = 3 * n - 5 + i
        col_start = 4 * i - 4
        A[row_idx][col_start] = 6 * X[i]
        A[row_idx][col_start + 1] = 2
        A[row_idx][col_start + 2] = 0
        A[row_idx][col_start + 3] = 0
        A[row_idx][col_start + 4] = -6 * X[i]
        A[row_idx][col_start + 5] = -2
        A[row_idx][col_start + 6] = 0
        A[row_idx][col_start + 7] = 0
        b[row_idx] = 0

    A[m - 2][0] = 6 * X[0]
    A[m - 2][1] = 2
    b[m - 2] = 0

    A[m - 1][m - 4] = 6 * X[-1]
    A[m - 1][m - 3] = 2
    b[m - 1] = 0

    Saux = partial_pivoting(A, b, return_stages=False)

    Coef = []
    for i in range(n - 1):
        Coef.append([Saux[4 * i], Saux[4 * i + 1], Saux[4 * i + 2], Saux[4 * i + 3]])

    return {
        "coefficients": Coef,
        "x_points": X,
        "y_points": Y,
        "num_segments": n - 1
    }

def format_matrix(matrix):
    return '\n'.join([' ' + '  '.join([f'{x:.6f}' for x in row]) + ' ' for row in matrix])

def format_factorization_result(method_name, stages, solution, has_P=False, P=None):
    result = f"{method_name}\n\nResults:\n"

    for idx, stage in enumerate(stages):
        result += f"\nStage {idx}\n\n"
        # Show augmented matrix
        result += format_matrix(stage['M']) + "\n"

        # Show L and U after first stage
        if idx > 0:
            result += "\nL:\n"
            result += format_matrix(stage['L']) + "\n"
            result += "\nU:\n"
            result += format_matrix(stage['U']) + "\n"
            if has_P and 'P' in stage:
                result += "\nP:\n"
                result += format_matrix(stage['P']) + "\n"

    result += "\n\nAfter applying forward and backward substitution\n\nx:\n"
    for x in solution:
        result += f"{x:.6f}\n"

    return result

@app.route('/')
def index():
    return render_template('index.html',
                         result=None,
                         iterations=None,
                         headers=None,
                         error=None,
                         method='bisection',
                         form_data={})

@app.route('/calculate', methods=['POST'])
def calculate():
    try:
        method = request.form['method']
        result = None
        iterations = []
        headers = []

        if method in ['bisection', 'false_position', 'incremental_search', 'fixed_point', 'newton', 'secant', 'newton_multiple']:
            # Methods for finding roots
            expr = request.form['function']
            f = create_function(expr)

            if method == 'incremental_search':
                x0 = float(request.form['x0'])
                delta = float(request.form['delta'])
                nmax = int(request.form['nmax'])
                res = incremental_search(f, x0, delta, nmax)

                result = "Incremental Search\n\nResults:\n\n"

                if res['intervals']:
                    for interval in res['intervals']:
                        a, b = interval
                        result += f"There is a root of f in [{a:.10f},{b:.10f}]\n"
                    result += f"_______________________________________________________"
                else:
                    result = "No roots found in the specified interval"

                iterations = []
                headers = []

            elif method == 'bisection':
                a, b = float(request.form['a']), float(request.form['b'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = bisection_method(f, a, b, tol, nmax)

                result = "Bisection\n\nResults table:\n\n"
                result += "| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\n"

                prev_xm = None
                for iter_data in res['iterations']:
                    i, a_val, b_val, c_val, fc_val = iter_data
                    error_str = ""
                    if prev_xm is not None:
                        error = abs(c_val - prev_xm)
                        error_str = f"|  {error:.1e}  |"
                    else:
                        error_str = "|           |"

                    result += f"|  {i+1:<2} | {a_val:.10f}  | {c_val:.10f}  |  {b_val:.10f}  |  {fc_val:.1e}  {error_str} \n"
                    prev_xm = c_val

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

            elif method == 'false_position':
                a, b = float(request.form['a']), float(request.form['b'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = false_position_method(f, a, b, tol, nmax)

                result = "False Position\n\nResults table:\n\n"
                result += "| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\n"

                prev_xm = None
                for iter_data in res['iterations']:
                    i, a_val, b_val, c_val, fc_val = iter_data
                    error_str = ""
                    if prev_xm is not None:
                        error = abs(c_val - prev_xm)
                        error_str = f"|  {error:.1e}  |"
                    else:
                        error_str = "|           |"

                    result += f"|  {i+1:<2} | {a_val:.10f}  | {c_val:.10f}  |  {b_val:.10f}  |  {fc_val:.1e}  {error_str} \n"
                    prev_xm = c_val

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

            elif method == 'fixed_point':
                g_expr = request.form['g_function']
                g = create_function(g_expr)
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = fixed_point_method(g, x0, tol, nmax)

                result = "Fixed Point\n\nResults table:\n\n"
                result += "| iter|     xi         |     g(xi)      |   f(xi)   |     E     |\n"

                prev_xi = x0
                for iter_data in res['iterations']:
                    i, x_old, x_new = iter_data
                    fx = f(x_new)
                    error = abs(x_new - x_old)
                    error_str = f"|  {error:.1e}  |" if i > 0 else "|           |"

                    result += f"|  {i:<2} | {x_old:.10f}  | {x_new:.10f}  |  {fx:.1e}  {error_str} \n"

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

            elif method == 'newton':
                df = create_function(request.form['derivative']) if request.form.get('derivative') else None
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = newton_method(f, df, x0, tol, nmax)

                result = "Newton\n\nResults table:\n\n"
                result += "| iter|     xi         |   f(xi)   |     E     |\n"

                result += f"|  0  | {x0:.10f}   |  {f(x0):.1e}  |           |  \n"
                for iter_data in res['iterations']:
                    i, x_val, fx_val, dfx_val, x_next = iter_data
                    error = abs(x_next - x_val)
                    result += f"|  {i+1:<2} | {x_next:.10f}  |  {f(x_next):.1e}  |  {error:.1e}  | \n"

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

            elif method == 'secant':
                x0, x1 = float(request.form['x0']), float(request.form['x1'])
                tol, nmax = float(request.form['tolerance']), int(request.form['nmax'])
                res = secant_method(f, x0, x1, tol, nmax)

                result = "Secant\n\nResults table:\n\n"
                result += "| iter|     xi         |   f(xi)   |     E     |\n"

                result += f"|  0  | {x0:.10f}   |  {f(x0):.1e}  |           |  \n"
                result += f"|  1  | {x1:.10f}   |  {f(x1):.1e}  |           |  \n"

                prev_x = x1
                for idx, iter_data in enumerate(res['iterations']):
                    i, x0_val, x1_val, x2_val, fx1_val = iter_data
                    error = abs(x2_val - prev_x)
                    result += f"|  {idx+2:<2} | {x2_val:.10f}  |  {f(x2_val):.1e}  |  {error:.1e}  | \n"
                    prev_x = x2_val

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

            elif method == 'newton_multiple':
                df = create_function(request.form['derivative']) if request.form.get('derivative') else None
                d2f = create_function(request.form['second_derivative']) if request.form.get('second_derivative') else None
                x0, tol, nmax = float(request.form['x0']), float(request.form['tolerance']), int(request.form['nmax'])
                res = newton_multiple_roots(f, df, d2f, x0, tol, nmax)

                result = "Multiple Roots\n\nResults table:\n\n"
                result += "| iter|     xi         |   f(xi)   |     E     |\n"

                result += f"|  0  | {x0:.10f}   |  {f(x0):.1e}  |           |  \n"
                for iter_data in res['iterations']:
                    i, x_val, fx_val, dfx_val, d2fx_val, x_next = iter_data
                    error = abs(x_next - x_val)
                    result += f"|  {i+1:<2} | {x_next:.10f}   |  {f(x_next):.1e}  |  {error:.1e}  | \n"

                result += f"\nFound an approximation of the root at {res['root']:.15f}"
                iterations = []
                headers = []

        elif method in ['gauss', 'partial_pivot', 'total_pivot', 'lu_simple', 'lu_partial', 'crout', 'doolittle', 'cholesky', 'jacobi', 'gauss_seidel', 'sor', 'forward_substitution', 'backward_substitution']:
            # Methods for linear systems
            if method in ['forward_substitution', 'backward_substitution']:
                # These methods receive augmented matrix [A|b]
                M = parse_matrix(request.form['augmented_matrix'])
            else:
                A = parse_matrix(request.form['matrix_A'])
                b = [float(x) for x in request.form['vector_b'].split()]

            if method in ['jacobi', 'gauss_seidel', 'sor']:
                # Iterative methods need x0, tol, nmax
                x0_str = request.form.get('x0_vector', '0 ' * len(b)).strip()
                x0 = [float(x) for x in x0_str.split()] if x0_str else [0.0] * len(b)
                tol = float(request.form.get('tolerance', 1e-6))
                nmax = int(request.form.get('nmax', 50))

                if method == 'jacobi':
                    res = jacobi(A, b, x0, tol, nmax)
                    method_name = "Jacobi"
                elif method == 'gauss_seidel':
                    res = gauss_seidel(A, b, x0, tol, nmax)
                    method_name = "Gauss-Seidel"
                else:  # sor
                    omega = float(request.form.get('omega', 1.5))
                    res = sor(A, b, x0, omega, tol, nmax)
                    method_name = f"SOR (relajación)"

                solution = res['solution']
                T = res['T']
                C = res['C']
                spectral_radius = res['spectral_radius']

                # Format output
                result = f"{method_name}\n\nResults:\n\nT:\n"
                result += format_matrix(T) + "\n"
                result += "\nC:\n"
                for c_val in C:
                    result += f" {c_val:.6f} "
                result += "\n\nSpectral radius:\n"
                result += f" {spectral_radius:.6f}\n\n"

                # Format iteration table
                result += "| iter|      E     |"
                for i in range(len(solution)):
                    result += " "
                result += "\n"

                result += f"|  0  |            | " + "  ".join([f"{x0[i]:.6f}" for i in range(len(solution))]) + " \n"

                for iter_data in res['iterations']:
                    k = iter_data[0]
                    error = iter_data[1]
                    x_vals = iter_data[2:]
                    result += f"|  {k+1:<2} |   {error:.1e}  | " + "  ".join([f"{x:.6f}" for x in x_vals]) + " \n"

                iterations = []
                headers = []

            elif method == 'gauss':
                res = gauss_elimination(A, b, return_stages=True)
                solution = res['solution']
                stages = res['stages']

                # Format stages
                result = "Simple Gaussian Elimination\n\nResults:\n"
                for i, (desc, matrix) in enumerate(stages):
                    result += f"\nStage {i}\n\n"
                    for row in matrix:
                        result += " " + "  ".join([f"{x:.6f}" for x in row]) + " \n"

                result += "\n\nAfter applying backward substitution\n\nx:\n"
                for x in solution:
                    result += f"{x:.6f}\n"
            elif method == 'partial_pivot':
                res = partial_pivoting(A, b, return_stages=True)
                solution = res['solution']
                stages = res['stages']

                # Format stages
                result = "Gaussian Elimination with Partial Pivoting\n\nResults:\n"
                for i, (desc, matrix) in enumerate(stages):
                    result += f"\nStage {i}\n\n"
                    for row in matrix:
                        result += " " + "  ".join([f"{x:.6f}" for x in row]) + " \n"

                result += "\n\nAfter applying backward substitution\n\nx:\n"
                for x in solution:
                    result += f"{x:.6f}\n"
            elif method == 'total_pivot':
                res = total_pivoting(A, b, return_stages=True)
                solution = res['solution']
                stages = res['stages']
                marks = res['marks']

                # Format stages
                result = "Gaussian Elimination with Total Pivoting\n\nResults:\n"
                for i, (desc, matrix) in enumerate(stages):
                    result += f"\nStage {i}\n\n"
                    for row in matrix:
                        result += " " + "  ".join([f"{x:.6f}" for x in row]) + " \n"

                result += "\n\nAfter applying backward substitution\n\nx:\n"
                for x in solution:
                    result += f"{x:.6f}\n"
            elif method == 'lu_simple':
                res = lu_simple(A, b)
                result = format_factorization_result("LU with Simple Gaussian:", res['stages'], res['solution'])
            elif method == 'lu_partial':
                res = lu_partial_pivot(A, b)
                # Use helper function but add P matrix
                result = "LU with Partial Pivoting\n\nResults:\n"

                for idx, stage in enumerate(res['stages']):
                    result += f"\nStage {idx}\n\n"
                    result += format_matrix(stage['M']) + "\n"

                    if idx > 0:
                        result += "\nL:\n"
                        result += format_matrix(stage['L']) + "\n"
                        result += "\nU:\n"
                        result += format_matrix(stage['U']) + "\n"
                        result += "\nP:\n"
                        result += format_matrix(stage['P']) + "\n"

                result += "\n\nAfter applying forward and backward substitution\n\nx:\n"
                for x in res['solution']:
                    result += f"{x:.6f}\n"
            elif method == 'crout':
                res = crout(A, b)
                result = format_factorization_result("Crout", res['stages'], res['solution'])
            elif method == 'doolittle':
                res = doolittle(A, b)
                result = format_factorization_result("Doolittle", res['stages'], res['solution'])
            elif method == 'cholesky':
                res = cholesky(A, b)
                result = format_factorization_result("Cholesky", res['stages'], res['solution'])
            elif method == 'forward_substitution':
                res = forward_substitution_method(M)
                solution = res['solution']
                stages = res['stages']

                # Format output
                result = "Forward Substitution\n\nResults:\n\n"

                # Show initial matrix
                desc, matrix = stages[0]
                result += f"{desc}:\n\n"
                for row in matrix:
                    result += " " + "  ".join([f"{x:.6f}" for x in row]) + " \n"

                # Show calculation steps
                result += "\nCalculation steps:\n\n"
                desc, steps = stages[1]
                for step in steps:
                    result += step + "\n"

                result += "\nSolution:\n\nx:\n"
                for x_val in solution:
                    result += f"{x_val:.6f}\n"
            elif method == 'backward_substitution':
                res = backward_substitution_method(M)
                solution = res['solution']
                stages = res['stages']

                # Format output
                result = "Backward Substitution\n\nResults:\n\n"

                # Show initial matrix
                desc, matrix = stages[0]
                result += f"{desc}:\n\n"
                for row in matrix:
                    result += " " + "  ".join([f"{x:.6f}" for x in row]) + " \n"

                # Show calculation steps
                result += "\nCalculation steps:\n\n"
                desc, steps = stages[1]
                for step in steps:
                    result += step + "\n"

                result += "\nSolution:\n\nx:\n"
                for x_val in solution:
                    result += f"{x_val:.6f}\n"

            if method not in ['gauss', 'partial_pivot', 'total_pivot', 'lu_simple', 'lu_partial', 'crout', 'doolittle', 'cholesky', 'jacobi', 'gauss_seidel', 'sor', 'forward_substitution', 'backward_substitution']:
                result = f"Solution: {', '.join([f'{x:.8g}' for x in solution])}"

        elif method in ['vandermonde', 'divided_differences', 'lagrange', 'linear_spline', 'quadratic_spline', 'cubic_spline']:
            # Interpolation methods
            X_str = request.form['x_values']
            Y_str = request.form['y_values']
            X = [float(x) for x in X_str.split()]
            Y = [float(y) for y in Y_str.split()]

            if len(X) != len(Y):
                raise ValueError("X and Y must have the same number of points")

            if method == 'vandermonde':
                res = vandermonde(X, Y)
                coef = res['coefficients']
                degree = res['degree']

                # Build Vandermonde matrix
                n = len(X)
                V = []
                for i in range(n):
                    row = [X[i]**(n - 1 - j) for j in range(n)]
                    V.append(row)

                # Format output
                result = "Vandermonde\n\nResults:\n\nVandermonde matrix:\n\n"
                result += format_matrix(V) + "\n"
                result += "\nPolynomial coefficients:\n\n"
                for c in coef:
                    result += f" {c:.6f} "
                result += "\n\nPolynomial:\n\n"

                # Format polynomial
                poly_terms = []
                for i in range(len(coef)):
                    exp = degree - i
                    if exp > 1:
                        poly_terms.append(f"{coef[i]:.6f}x^{exp}")
                    elif exp == 1:
                        poly_terms.append(f"{coef[i]:.6f}x")
                    else:
                        poly_terms.append(f"{coef[i]:.6f}")

                result += "+".join(poly_terms).replace("+-", "-")

            elif method == 'divided_differences':
                res = divided_differences(X, Y)
                coef = res['coefficients']
                degree = res['degree']
                x_points = res['x_points']
                table = res['table']

                # Format output
                result = "Newton\n\nResults:\n\nDivided differences table:\n\n"
                result += format_matrix(table) + "\n"
                result += "\nNewton polynomial coefficients:\n\n"
                for c in coef:
                    result += f" {c:.6f} "
                result += "\n\nNewton polynomial:\n\n"

                # Format Newton polynomial
                poly_parts = [f"{coef[0]:.6f}"]
                for i in range(1, len(coef)):
                    # Build factors as (x - x_j), not (x + sign*x_j)
                    factors = "".join([f"(x-({x_points[j]:.6f}))" if x_points[j] != 0 else "x" for j in range(i)])
                    poly_parts.append(f"{coef[i]:.6f}{factors}")

                # Simplify: (x-(-a)) becomes (x+a), (x-(0)) becomes x
                result += "+".join(poly_parts).replace("+-", "-").replace("(x-(-", "(x+").replace("))", ")")

            elif method == 'lagrange':
                res = lagrange(X, Y)
                coef = res['coefficients']
                degree = res['degree']
                L_basis = res['lagrange_basis']

                # Format output
                result = "Lagrange\n\nResults:\n\nLagrange interpolating polynomials:\n\n"

                # Format Lagrange basis polynomials
                for i in range(len(L_basis)):
                    Li = L_basis[i]
                    Li_degree = len(Li) - 1
                    terms = []
                    for j in range(len(Li)):
                        exp = Li_degree - j
                        # Show all terms, even near-zero ones
                        if exp > 1:
                            terms.append(f"{Li[j]:.6f}x^{exp}")
                        elif exp == 1:
                            terms.append(f"{Li[j]:.6f}x")
                        else:
                            terms.append(f"{Li[j]:.6f}")
                    result += "+".join(terms).replace("+-", "-") + f"   //L{i}\n"

                result += "\n\nPolynomial:\n\n"
                # Show combination
                poly_parts = [f"{Y[i]:.1f}*L{i}" for i in range(len(Y))]
                result += "+".join(poly_parts).replace("+-", "-")

            elif method == 'linear_spline':
                res = linear_spline(X, Y)
                coef = res['coefficients']
                num_segments = res['num_segments']
                x_points = res['x_points']

                # Format output
                result = "Linear Splines\n\nResults:\n\nSpline coefficients:\n\n"
                for i in range(num_segments):
                    a_i, b_i = coef[i]
                    result += f" {a_i:.6f}  {b_i:.6f} \n"

                result += "\nSplines:\n\n"
                for i in range(num_segments):
                    a_i, b_i = coef[i]
                    sign = "+" if b_i >= 0 else ""
                    result += f"{a_i:.6f}x{sign}{b_i:.6f}\n"

            elif method == 'quadratic_spline':
                res = quadratic_spline(X, Y)
                coef = res['coefficients']
                num_segments = res['num_segments']
                x_points = res['x_points']

                # Format output
                result = "Quadratic Splines\n\nResults:\n\nSpline coefficients:\n\n"
                for i in range(num_segments):
                    a_i, b_i, c_i = coef[i]
                    result += f" {a_i:.6f}  {b_i:.6f}  {c_i:.6f} \n"

                result += "\nSplines:\n\n"
                for i in range(num_segments):
                    a_i, b_i, c_i = coef[i]
                    b_sign = "+" if b_i >= 0 else ""
                    c_sign = "+" if c_i >= 0 else ""
                    result += f"{a_i:.6f}x^2{b_sign}{b_i:.6f}x{c_sign}{c_i:.6f}\n"

            elif method == 'cubic_spline':
                res = cubic_spline(X, Y)
                coef = res['coefficients']
                num_segments = res['num_segments']
                x_points = res['x_points']

                # Format output
                result = "Cubic Splines\n\nResults:\n\nSpline coefficients:\n\n"
                for i in range(num_segments):
                    a_i, b_i, c_i, d_i = coef[i]
                    result += f" {a_i:.6f}  {b_i:.6f}  {c_i:.6f}  {d_i:.6f} \n"

                result += "\nSplines:\n\n"
                for i in range(num_segments):
                    a_i, b_i, c_i, d_i = coef[i]
                    b_sign = "+" if b_i >= 0 else ""
                    c_sign = "+" if c_i >= 0 else ""
                    d_sign = "+" if d_i >= 0 else ""
                    result += f"{a_i:.6f}x^3{b_sign}{b_i:.6f}x^2{c_sign}{c_i:.6f}x{d_sign}{d_i:.6f}\n"

        return render_template('index.html',
                             result=result,
                             iterations=iterations,
                             headers=headers,
                             error=None,
                             method=method,
                             form_data=request.form)

    except Exception as e:
        return render_template('index.html',
                             result=None,
                             iterations=None,
                             headers=None,
                             error=str(e),
                             method=request.form.get('method', 'bisection'),
                             form_data=request.form)

if __name__ == '__main__':
    app.run(debug=True)