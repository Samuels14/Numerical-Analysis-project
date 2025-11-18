# Numerical Calculator - JavaScript Version

Calculadora de análisis numérico implementada en JavaScript con Node.js y Express.

## Instalación

```bash
npm install
```

## Ejecución

```bash
npm start
```

O para desarrollo con auto-reload:

```bash
npm run dev
```

La aplicación estará disponible en: http://localhost:3000

## Librerías Utilizadas

- **Express**: Framework web para Node.js
- **mathjs**: Librería para evaluación de expresiones matemáticas
- **numeric**: Librería para cálculos numéricos (matrices, álgebra lineal)
- **EJS**: Motor de plantillas para renderizar HTML

## Métodos Implementados

### Ecuaciones No Lineales
- Búsqueda Incremental
- Bisección
- Falsa Posición
- Punto Fijo
- Newton-Raphson
- Secante
- Newton Raíces Múltiples

### Sistemas Lineales
- Eliminación Gaussiana (Simple, Pivoteo Parcial, Pivoteo Total)
- Factorización LU (Simple, Pivoteo Parcial)
- Crout
- Doolittle
- Cholesky
- Jacobi
- Gauss-Seidel
- SOR

### Interpolación
- Vandermonde
- Diferencias Divididas
- Lagrange
- Splines (Lineal, Cuadrático, Cúbico)

### Métodos Adicionales
- Sustitución Progresiva
- Sustitución Regresiva
