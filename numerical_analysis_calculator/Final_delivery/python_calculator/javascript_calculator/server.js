const express = require('express');
const path = require('path');
const methods = require('./lib/methods');

const app = express();
const PORT = process.env.PORT || 3000;

// Middleware
app.use(express.urlencoded({ extended: true }));
app.use(express.json());
app.use('/static', express.static(path.join(__dirname, 'static')));

// View engine
app.set('view engine', 'ejs');
app.set('views', path.join(__dirname, 'views'));

// Routes
app.get('/', (req, res) => {
    res.render('index', {
        result: null,
        iterations: null,
        headers: null,
        error: null,
        method: 'bisection',
        form_data: {},
        warnings: [],
        interpolation_data: null
    });
});

app.post('/calculate', async (req, res) => {
    try {
        const method = req.body.method;
        let result = null;
        let iterations = [];
        let headers = [];
        let warnings = [];
        let interpolation_data = null;

        // Call the appropriate method
        const methodResult = await methods.calculate(method, req.body);
        
        result = methodResult.result;
        iterations = methodResult.iterations || [];
        headers = methodResult.headers || [];
        warnings = methodResult.warnings || [];

        // Prepare interpolation data for graphing
        if (['vandermonde', 'divided_differences', 'lagrange', 'linear_spline', 'quadratic_spline', 'cubic_spline'].includes(method)) {
            const xPoints = req.body.x_values.trim().split(/\s+/).map(parseFloat);
            const yPoints = req.body.y_values.trim().split(/\s+/).map(parseFloat);
            
            interpolation_data = {
                x_points: xPoints,
                y_points: yPoints,
                method: method
            };
        }

        res.render('index', {
            result,
            iterations,
            headers,
            error: null,
            method,
            form_data: req.body,
            warnings,
            interpolation_data
        });

    } catch (error) {
        res.render('index', {
            result: null,
            iterations: null,
            headers: null,
            error: error.message,
            method: req.body.method || 'bisection',
            form_data: req.body,
            warnings: [],
            interpolation_data: null
        });
    }
});

app.listen(PORT, () => {
    console.log('╔════════════════════════════════════════════════╗');
    console.log('║  Numerical Calculator - JavaScript Version    ║');
    console.log('╚════════════════════════════════════════════════╝');
    console.log('');
    console.log(`✅ Server running on: http://localhost:${PORT}`);
    console.log('✅ Press Ctrl+C to stop');
    console.log('');
});
