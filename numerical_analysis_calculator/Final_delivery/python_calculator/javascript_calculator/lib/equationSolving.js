const { createFunction, safeEval, numericDerivative, numericSecondDerivative } = require('./utils');

// Incremental Search
function incremental_search(formData) {
    const f = createFunction(formData.function);
    const x0 = parseFloat(formData.x0);
    const delta = parseFloat(formData.delta);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const intervals = [];
    const suggestions = [];
    
    if (nmax <= 0) {
        throw new Error("Nmax must be a positive integer.");
    }
    if (delta === 0) {
        throw new Error("Increment Δ cannot be zero. Use a positive or negative non-zero value.");
    }
    
    let x_prev = x0;
    let f_prev = safeEval(f, x_prev);
    let x_curr = x_prev + delta;
    let f_curr = safeEval(f, x_curr);
    
    // Validate first evaluations
    if (!isFinite(f_prev) || !isFinite(f_curr)) {
        throw new Error(
            "Initial evaluation of f(x) yielded NaN or Inf. " +
            "Check function definition or adjust the starting interval."
        );
    }
    
    for (let i = 0; i < nmax; i++) {
        iterations.push([i, x_prev, f_prev, x_curr, f_curr]);
        
        if (f_prev * f_curr < 0) {
            intervals.push([x_prev, x_curr]);
        }
        
        x_prev = x_curr;
        f_prev = f_curr;
        x_curr = x_prev + delta;
        f_curr = safeEval(f, x_curr);
        
        if (isNaN(f_curr) || !isFinite(f_curr)) {
            suggestions.push(
                `At iteration ${i+1}, f(x) produced a non-finite value. ` +
                "The function may not be well-defined over the explored domain."
            );
            break;
        }
    }
    
    let result = "Incremental Search\\n\\nResults:\\n\\n";
    if (intervals.length > 0) {
        for (const interval of intervals) {
            result += `There is at least one root of f in [${interval[0].toFixed(10)},${interval[1].toFixed(10)}]\n`;
        }
        result += "_______________________________________________________";
    } else {
        result = "No sign-change intervals were found.\n";
        suggestions.push(
            "No sign changes were found within the explored range. " +
            "Consider increasing Nmax, adjusting Δ, or checking if the function actually has roots nearby."
        );
    }
    
    return { result, iterations: [], headers: [], warnings: suggestions };
}

// Bisection Method
function bisection(formData) {
    const f = createFunction(formData.function);
    let a = parseFloat(formData.a);
    let b = parseFloat(formData.b);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const warnings = [];
    
    if (a >= b) {
        throw new Error("Invalid interval: bisection requires a < b.");
    }
    if (tol <= 0) {
        throw new Error("Tolerance must be a positive number.");
    }
    if (nmax <= 0) {
        throw new Error("Nmax must be a positive integer.");
    }
    
    let fa = safeEval(f, a);
    let fb = safeEval(f, b);
    
    if (!isFinite(fa) || !isFinite(fb)) {
        throw new Error(
            "f(a) or f(b) resulted in NaN or Inf. " +
            "Check the function or the interval bounds."
        );
    }
    
    if (fa * fb > 0) {
        throw new Error(
            "Function does not change sign over [a, b]. " +
            "Bisection requires f(a) * f(b) < 0."
        );
    }
    
    let pm = (a + b) / 2;
    let fpm = safeEval(f, pm);
    let E = 1000;
    let cnt = 1;
    
    let result = "Bisection\\n\\nResults table:\\n\\n";
    result += "| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\\n";
    
    let prev_pm = null;
    iterations.push([0, a, b, pm, fpm]);
    
    while (E > tol && cnt < nmax) {
        if (fa * fpm < 0) {
            b = pm;
            fb = fpm;
        } else {
            a = pm;
            fa = fpm;
        }
        
        prev_pm = pm;
        pm = (a + b) / 2;
        fpm = safeEval(f, pm);
        E = Math.abs(pm - prev_pm);
        
        iterations.push([cnt, a, b, pm, fpm]);
        cnt++;
    }
    
    for (const iter of iterations) {
        const [i, a_val, b_val, c_val, fc_val] = iter;
        const errorStr = i > 0 ? `|  ${Math.abs(c_val - iterations[i-1][3]).toExponential(1)}  |` : "|           |";
        result += `|  ${(i+1).toString().padEnd(2)} | ${a_val.toFixed(10)}  | ${c_val.toFixed(10)}  |  ${b_val.toFixed(10)}  |  ${fc_val.toExponential(1)}  ${errorStr} \\n`;
    }
    
    result += `\nFound an approximation of the root at ${pm.toFixed(15)}`;
    
    const converged = E <= tol;
    if (!converged) {
        warnings.push(
            "Bisection did not reach the desired tolerance within Nmax iterations. " +
            "Try increasing Nmax or using a larger tolerance."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// False Position Method
function false_position(formData) {
    const f = createFunction(formData.function);
    let a = parseFloat(formData.a);
    let b = parseFloat(formData.b);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const warnings = [];
    
    if (a >= b) {
        throw new Error("Invalid interval: false position requires a < b.");
    }
    if (tol <= 0) {
        throw new Error("Tolerance must be a positive number.");
    }
    if (nmax <= 0) {
        throw new Error("Nmax must be a positive integer.");
    }
    
    let fa = safeEval(f, a);
    let fb = safeEval(f, b);
    
    if (!isFinite(fa) || !isFinite(fb)) {
        throw new Error(
            "f(a) or f(b) resulted in NaN or Inf. " +
            "Check the function or the interval bounds."
        );
    }
    
    if (fa * fb > 0) {
        throw new Error(
            "Function does not change sign over [a, b]. " +
            "False position requires f(a) * f(b) < 0."
        );
    }
    
    if (fb - fa === 0) {
        throw new Error(
            "Cannot apply false position because f(a) and f(b) are equal. " +
            "Try modifying the interval or using bisection."
        );
    }
    
    let pm = (fb * a - fa * b) / (fb - fa);
    let fpm = safeEval(f, pm);
    let E = 1000;
    let cnt = 1;
    
    const iterations = [];
    iterations.push([0, a, b, pm, fpm]);
    
    while (E > tol && cnt < nmax) {
        if (fa * fpm < 0) {
            b = pm;
            fb = fpm;
        } else {
            a = pm;
            fa = fpm;
        }
        
        const prev_pm = pm;
        fa = safeEval(f, a);
        fb = safeEval(f, b);
        
        if (fb - fa === 0) {
            warnings.push(
                "During iterations, f(a) ≈ f(b), causing instability in the false position formula. " +
                "Consider switching to bisection for this interval."
            );
            break;
        }
        
        pm = (fb * a - fa * b) / (fb - fa);
        fpm = safeEval(f, pm);
        E = Math.abs(pm - prev_pm);
        
        iterations.push([cnt, a, b, pm, fpm]);
        cnt++;
    }
    
    let result = "False Position\\n\\nResults table:\\n\\n";
    result += "| iter|     a         |     xm        |     b          |     f(Xm)  |     E     |\\n";
    
    for (const iter of iterations) {
        const [i, a_val, b_val, c_val, fc_val] = iter;
        const errorStr = i > 0 ? `|  ${Math.abs(c_val - iterations[i-1][3]).toExponential(1)}  |` : "|           |";
        result += `|  ${(i+1).toString().padEnd(2)} | ${a_val.toFixed(10)}  | ${c_val.toFixed(10)}  |  ${b_val.toFixed(10)}  |  ${fc_val.toExponential(1)}  ${errorStr} \\n`;
    }
    
    result += `\nFound an approximation of the root at ${pm.toFixed(15)}`;
    
    const converged = E <= tol;
    if (!converged) {
        warnings.push(
            "False position did not reach the desired tolerance within Nmax iterations. " +
            "This method often stalls; consider switching to bisection or adjusting the interval."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Fixed Point Method
function fixed_point(formData) {
    const f = createFunction(formData.function);
    const g = createFunction(formData.g_function);
    let x0 = parseFloat(formData.x0);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const warnings = [];
    
    if (tol <= 0) {
        throw new Error("Tolerance must be a positive number.");
    }
    if (nmax <= 0) {
        throw new Error("Nmax must be a positive integer.");
    }
    let converged = false;
    
    for (let i = 0; i < nmax; i++) {
        const x1 = safeEval(g, x0, "g");
        
        if (isNaN(x1) || !isFinite(x1)) {
            warnings.push(
                `At iteration ${i}, g(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "Check the g function definition and the current iterate."
            );
            break;
        }
        
        iterations.push([i, x0, x1]);
        
        if (Math.abs(x1 - x0) < tol) {
            converged = true;
            x0 = x1;
            break;
        }
        
        x0 = x1;
    }
    
    let result = "Fixed Point\\n\\nResults table:\\n\\n";
    result += "| iter|     xi         |     g(xi)      |   f(xi)   |     E     |\\n";
    
    for (const iter of iterations) {
        const [i, x_old, x_new] = iter;
        const fx = safeEval(f, x_new);
        const error = Math.abs(x_new - x_old);
        const errorStr = i > 0 ? `|  ${error.toExponential(1)}  |` : "|           |";
        result += `|  ${i.toString().padEnd(2)} | ${x_old.toFixed(10)}  | ${x_new.toFixed(10)}  |  ${fx.toExponential(1)}  ${errorStr} \\n`;
    }
    
    result += `\nFound an approximation of the root at ${x0.toFixed(15)}`;
    
    if (!converged) {
        warnings.push(
            "Fixed point iteration did not reach the desired tolerance within Nmax iterations. " +
            "Try a different g(x) function, a different initial guess, or increase Nmax. " +
            "Ensure |g'(x)| < 1 near the root for convergence."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Newton's Method
function newton(formData) {
    const f = createFunction(formData.function);
    const df = formData.derivative ? createFunction(formData.derivative) : null;
    let x0 = parseFloat(formData.x0);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const warnings = [];
    let converged = false;
    
    for (let i = 0; i < nmax; i++) {
        const fx = safeEval(f, x0);
        
        if (isNaN(fx) || !isFinite(fx)) {
            warnings.push(
                `At iteration ${i}, f(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "Check the function definition and the current iterate."
            );
            break;
        }
        
        const dfx = df ? df(x0) : numericDerivative(f, x0);
        
        if (isNaN(dfx) || !isFinite(dfx)) {
            warnings.push(
                `At iteration ${i}, f'(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "The derivative may not be well-defined at this point."
            );
            break;
        }
        
        if (Math.abs(dfx) < 1e-14) {
            warnings.push(
                `At iteration ${i}, the derivative f'(x) ≈ 0 at x = ${x0}. ` +
                "Newton's method cannot proceed reliably when the derivative is (almost) zero. " +
                "Try a different initial guess or use a bracketing method such as bisection."
            );
            break;
        }
        
        const x1 = x0 - fx / dfx;
        iterations.push([i, x0, fx, dfx, x1]);
        
        if (Math.abs(x1 - x0) < tol) {
            converged = true;
            x0 = x1;
            break;
        }
        
        x0 = x1;
    }
    
    let result = "Newton\\n\\nResults table:\\n\\n";
    result += "| iter|     xi         |   f(xi)   |     E     |\\n";
    result += `|  0  | ${parseFloat(formData.x0).toFixed(10)}   |  ${safeEval(f, parseFloat(formData.x0)).toExponential(1)}  |           |  \\n`;
    
    for (const iter of iterations) {
        const [i, x_val, fx_val, dfx_val, x_next] = iter;
        const error = Math.abs(x_next - x_val);
        result += `|  ${(i+1).toString().padEnd(2)} | ${x_next.toFixed(10)}  |  ${safeEval(f, x_next).toExponential(1)}  |  ${error.toExponential(1)}  | \\n`;
    }
    
    result += `\nFound an approximation of the root at ${x0.toFixed(15)}`;
    
    if (!converged) {
        warnings.push(
            "Newton's method did not reach the desired tolerance within Nmax iterations. " +
            "Consider changing the initial guess, providing an analytical derivative if you are using the numerical one, " +
            "or switching to a more robust method like bisection or the secant method."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Secant Method
function secant(formData) {
    const f = createFunction(formData.function);
    let x0 = parseFloat(formData.x0);
    let x1 = parseFloat(formData.x1);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const warnings = [];
    let converged = false;
    
    if (x0 === x1) {
        throw new Error("Initial points x0 and x1 must be different for the secant method.");
    }
    
    for (let i = 0; i < nmax; i++) {
        const fx0 = safeEval(f, x0);
        const fx1 = safeEval(f, x1);
        
        if (isNaN(fx0) || !isFinite(fx0) || isNaN(fx1) || !isFinite(fx1)) {
            warnings.push(
                `At iteration ${i}, f(x) evaluated at x0 = ${x0} or x1 = ${x1} produced NaN or Inf. ` +
                "Check the function and the chosen initial points."
            );
            break;
        }
        
        const denom = fx1 - fx0;
        if (Math.abs(denom) < 1e-14) {
            warnings.push(
                `At iteration ${i}, f(x1) - f(x0) ≈ 0, which makes the secant update numerically unstable. ` +
                "Try different initial points or use a bracketing method."
            );
            break;
        }
        
        const x2 = x1 - fx1 * (x1 - x0) / denom;
        iterations.push([i, x0, x1, x2, fx1]);
        
        if (Math.abs(x2 - x1) < tol) {
            converged = true;
            x1 = x2;
            break;
        }
        
        x0 = x1;
        x1 = x2;
    }
    
    let result = "Secant\\n\\nResults table:\\n\\n";
    result += "| iter|     xi         |   f(xi)   |     E     |\\n";
    result += `|  0  | ${parseFloat(formData.x0).toFixed(10)}   |  ${safeEval(f, parseFloat(formData.x0)).toExponential(1)}  |           |  \\n`;
    result += `|  1  | ${parseFloat(formData.x1).toFixed(10)}   |  ${safeEval(f, parseFloat(formData.x1)).toExponential(1)}  |           |  \\n`;
    
    let prev_x = parseFloat(formData.x1);
    for (const [idx, iter] of iterations.entries()) {
        const [i, x0_val, x1_val, x2_val, fx1_val] = iter;
        const error = Math.abs(x2_val - prev_x);
        result += `|  ${(idx+2).toString().padEnd(2)} | ${x2_val.toFixed(10)}  |  ${safeEval(f, x2_val).toExponential(1)}  |  ${error.toExponential(1)}  | \\n`;
        prev_x = x2_val;
    }
    
    result += `\nFound an approximation of the root at ${x1.toFixed(15)}`;
    
    if (!converged) {
        warnings.push(
            "The secant method did not reach the desired tolerance before Nmax. " +
            "Try different initial points or a bracketing method (like bisection) if you know a sign-change interval."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

// Newton Multiple Roots
function newton_multiple(formData) {
    const f = createFunction(formData.function);
    const df = formData.derivative ? createFunction(formData.derivative) : null;
    const d2f = formData.second_derivative ? createFunction(formData.second_derivative) : null;
    let x0 = parseFloat(formData.x0);
    const tol = parseFloat(formData.tolerance);
    const nmax = parseInt(formData.nmax);
    
    const iterations = [];
    const warnings = [];
    
    if (tol <= 0) {
        throw new Error("Tolerance must be a positive number.");
    }
    if (nmax <= 0) {
        throw new Error("Nmax must be a positive integer.");
    }
    let converged = false;
    
    for (let i = 0; i < nmax; i++) {
        const fx = safeEval(f, x0);
        
        if (isNaN(fx) || !isFinite(fx)) {
            warnings.push(
                `At iteration ${i}, f(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "Check the function definition and the current iterate."
            );
            break;
        }
        
        const dfx = df ? df(x0) : numericDerivative(f, x0);
        
        if (isNaN(dfx) || !isFinite(dfx)) {
            warnings.push(
                `At iteration ${i}, f'(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "The derivative may not be well-defined at this point."
            );
            break;
        }
        
        const d2fx = d2f ? d2f(x0) : numericSecondDerivative(f, x0);
        
        if (isNaN(d2fx) || !isFinite(d2fx)) {
            warnings.push(
                `At iteration ${i}, f''(x) evaluated at x = ${x0} produced NaN or Inf. ` +
                "The second derivative may not be well-defined at this point."
            );
            break;
        }
        
        const denom = dfx * dfx - fx * d2fx;
        if (Math.abs(denom) < 1e-14) {
            warnings.push(
                `At iteration ${i}, the denominator [f'(x)]² - f(x)f''(x) ≈ 0 at x = ${x0}. ` +
                "The method cannot proceed reliably. Try a different initial guess or method."
            );
            break;
        }
        
        const x1 = x0 - (fx * dfx) / denom;
        iterations.push([i, x0, fx, dfx, d2fx, x1]);
        
        if (Math.abs(x1 - x0) < tol) {
            converged = true;
            x0 = x1;
            break;
        }
        
        x0 = x1;
    }
    
    let result = "Multiple Roots\\n\\nResults table:\\n\\n";
    result += "| iter|     xi         |   f(xi)   |     E     |\\n";
    result += `|  0  | ${parseFloat(formData.x0).toFixed(10)}   |  ${safeEval(f, parseFloat(formData.x0)).toExponential(1)}  |           |  \\n`;
    
    for (const iter of iterations) {
        const [i, x_val, fx_val, dfx_val, d2fx_val, x_next] = iter;
        const error = Math.abs(x_next - x_val);
        result += `|  ${(i+1).toString().padEnd(2)} | ${x_next.toFixed(10)}   |  ${safeEval(f, x_next).toExponential(1)}  |  ${error.toExponential(1)}  | \\n`;
    }
    
    result += `\nFound an approximation of the root at ${x0.toFixed(15)}`;
    
    if (!converged) {
        warnings.push(
            "The multiple-roots Newton method did not reach the desired tolerance. " +
            "Check the function and the initial point. This method requires good initial guesses near multiple roots."
        );
    }
    
    return { result, iterations: [], headers: [], warnings };
}

module.exports = {
    incremental_search,
    bisection,
    false_position,
    fixed_point,
    newton,
    secant,
    newton_multiple
};
