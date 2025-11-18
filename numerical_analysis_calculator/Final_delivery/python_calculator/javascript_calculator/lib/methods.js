const equationSolving = require('./equationSolving');
const linearSystems = require('./linearSystems');
const interpolation = require('./interpolation');
const substitution = require('./substitution');

async function calculate(method, formData) {
    // Equation solving methods
    const equationMethods = [
        'incremental_search', 'bisection', 'false_position', 
        'fixed_point', 'newton', 'secant', 'newton_multiple'
    ];
    
    // Linear system methods
    const linearMethods = [
        'gauss', 'partial_pivot', 'total_pivot',
        'lu_simple', 'lu_partial', 'crout', 'doolittle', 'cholesky',
        'jacobi', 'gauss_seidel', 'sor'
    ];
    
    // Interpolation methods
    const interpolationMethods = [
        'vandermonde', 'divided_differences', 'lagrange',
        'linear_spline', 'quadratic_spline', 'cubic_spline'
    ];
    
    // Substitution methods
    const substitutionMethods = [
        'forward_substitution', 'backward_substitution'
    ];

    if (equationMethods.includes(method)) {
        return equationSolving[method](formData);
    } else if (linearMethods.includes(method)) {
        return linearSystems[method](formData);
    } else if (interpolationMethods.includes(method)) {
        return interpolation[method](formData);
    } else if (substitutionMethods.includes(method)) {
        return substitution[method](formData);
    } else {
        throw new Error(`Unknown method: ${method}`);
    }
}

module.exports = { calculate };
