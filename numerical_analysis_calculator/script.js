/* Numerical Calculator - JavaScript */

// UTILITIES
const $ = sel => document.querySelector(sel)
const methods = document.querySelectorAll('.method-btn')
let current = 'bisection'

function setActive(btn){
  methods.forEach(b=>b.classList.remove('active'))
  btn.classList.add('active')
  current = btn.dataset.method
  renderParams()
}

methods.forEach(b=>b.addEventListener('click',()=>setActive(b)))

function renderParams(){
  const cont = $('#params')
  cont.innerHTML = ''
  // common parameters for root-finding methods
  if(['bisection','false_position','secant','incremental_search','newton','newton_multiple','fixed_point'].includes(current)){
    if(current==='incremental_search'){
      cont.innerHTML = `
        <label>x0</label><input id='p_x0' type='number' value='0' step='any'>
        <label>delta</label><input id='p_delta' type='number' value='0.5' step='any'>
        <label>Nmax</label><input id='p_Nmax' type='number' value='50'>
      `
    }else if(current==='bisection' || current==='false_position'){
      cont.innerHTML = `
        <label>a</label><input id='p_a' type='number' value='1' step='any'>
        <label>b</label><input id='p_b' type='number' value='2' step='any'>
        <div class='two-col'>
          <div><label>tolerance</label><input id='p_tol' type='number' value='1e-6' step='any'></div>
          <div><label>Nmax</label><input id='p_Nmax' type='number' value='50'></div>
        </div>
      `
    }else if(current==='secant'){
      cont.innerHTML = `
        <label>x0</label><input id='p_x0' type='number' value='1' step='any'>
        <label>x1</label><input id='p_x1' type='number' value='2' step='any'>
        <div class='two-col'>
          <div><label>tolerance</label><input id='p_tol' type='number' value='1e-6' step='any'></div>
          <div><label>Nmax</label><input id='p_Nmax' type='number' value='50'></div>
        </div>
      `
    }else if(current==='newton' || current==='newton_multiple'){
      cont.innerHTML = `
        <label>x0</label><input id='p_x0' type='number' value='1.5' step='any'>
        <div class='two-col'>
          <div><label>tolerance</label><input id='p_tol' type='number' value='1e-6' step='any'></div>
          <div><label>Nmax</label><input id='p_Nmax' type='number' value='50'></div>
        </div>
        <label class='small'>If you want, you can enter g'(x) or f''(x) symbolically using math.js; otherwise numerical derivatives are used.</label>
        <label>Symbolic derivative (optional)</label>
        <input id='p_df' type='text' placeholder='e.g: 3*x^2 - 1, cos(x) - 1/2'>
        <label>Second derivative (only for multiple roots)</label>
        <input id='p_d2f' type='text' placeholder='e.g: 6*x'>
      `
    }else if(current==='fixed_point'){
      cont.innerHTML = `
        <label>g(x)</label><input id='p_g' type='text' placeholder='e.g: (x^3 - 2)/x, cos(x)'>
        <label>x0</label><input id='p_x0' type='number' value='1' step='any'>
        <div class='two-col'>
          <div><label>tolerance</label><input id='p_tol' type='number' value='1e-6' step='any'></div>
          <div><label>Nmax</label><input id='p_Nmax' type='number' value='50'></div>
        </div>
      `
    }
  }

  // parameters for linear systems
  if(['gauss','partial_pivot','total_pivot'].includes(current)){
    cont.innerHTML = `
      <label>Size n</label><input id='p_n' type='number' value='3'>
      <label>A (each row in a line, separate coefficients with spaces)</label>
      <textarea id='p_A' placeholder='e.g: 2 3 1\n1 -1 2\n3 0 -1'>2 3 -1\n1 -2 3\n0 1 1</textarea>
      <label>b (values separated by spaces)</label>
      <input id='p_b' type='text' value='1 -4 2'>
    `
  }
}

renderParams()

// PARSE function safely using math.js
function parseF(expr){
  try{
    // Replace common aliases to make functions more intuitive
    expr = expr.replace(/\bln\(/g, 'log(')  // ln(x) -> log(x)
    expr = expr.replace(/\bLog\(/g, 'log(')  // Log(x) -> log(x)
    
    const node = math.parse(expr)
    const code = node.compile()
    return x => {
      const scope = {x: x}
      return code.evaluate(scope)
    }
  }catch(e){
    throw new Error('Error parsing function: '+e.message)
  }
}

function numericDerivative(f, x, h=1e-6){
  return (f(x+h)-f(x-h))/(2*h)
}
function numericSecondDerivative(f,x,h=1e-5){
  return (f(x+h)-2*f(x)+f(x-h))/(h*h)
}

// NUMERICAL METHODS
function incremental_search(f,x0,delta,Nmax){
  const it=[]
  let fx0 = f(x0)
  for(let i=0;i<Nmax;i++){
    const x1 = x0+delta
    const fx1 = f(x1)
    it.push([i,x0,fx0,x1,fx1])
    if(fx0*fx1<0) return {interval:[x0,x1],it}
    x0= x1; fx0=fx1
  }
  return {interval:null,it}
}

function bisection(f,a,b,tol,Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const c=(a+b)/2; const fc=f(c)
    it.push([i,a,b,c,fc])
    if(Math.abs(fc)<tol || (b-a)/2<tol) return {root:c,it}
    if(f(a)*fc<0) b=c; else a=c
  }
  return {root:(a+b)/2,it}
}

function false_position(f,a,b,tol,Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const fa=f(a), fb=f(b)
    const c=(a*fb - b*fa)/(fb-fa)
    const fc=f(c)
    it.push([i,a,b,c,fc])
    if(Math.abs(fc)<tol) return {root:c,it}
    if(fa*fc<0) b=c; else a=c
  }
  return {root:c,it}
}

function fixed_point(g,x0,tol,Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const x1=g(x0)
    it.push([i,x0,x1])
    if(Math.abs(x1-x0)<tol) return {root:x1,it}
    x0=x1
  }
  return {root:x0,it}
}

function newton_method(f, df, x0, tol, Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const fx=f(x0)
    const dfx = df ? df(x0) : numericDerivative(f,x0)
    if(dfx===0) break
    const x1 = x0 - fx/dfx
    it.push([i,x0,fx,dfx,x1])
    if(Math.abs(x1-x0)<tol) return {root:x1,it}
    x0=x1
  }
  return {root:x0,it}
}

function secant(f,x0,x1,tol,Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const fx0=f(x0), fx1=f(x1)
    if(fx1-fx0===0) break
    const x2 = x1 - fx1*(x1-x0)/(fx1-fx0)
    it.push([i,x0,x1,x2,fx1])
    if(Math.abs(x2-x1)<tol) return {root:x2,it}
    x0=x1; x1=x2
  }
  return {root:x1,it}
}

function newton_multiple(f, df, d2f, x0, tol, Nmax){
  const it=[]
  for(let i=0;i<Nmax;i++){
    const fx=f(x0)
    const dfx = df ? df(x0) : numericDerivative(f,x0)
    const d2fx = d2f ? d2f(x0) : numericSecondDerivative(f,x0)
    const denom = dfx*dfx - fx*d2fx
    if(denom===0) break
    const x1 = x0 - (fx*dfx)/denom
    it.push([i,x0,fx,dfx,d2fx,x1])
    if(Math.abs(x1-x0)<tol) return {root:x1,it}
    x0=x1
  }
  return {root:x0,it}
}

// LINEAR SYSTEMS
function parseMatrix(Atext){
  const rows = Atext.trim().split('\n').map(r=>r.trim()).filter(r=>r.length)
  return rows.map(r=>r.split(/\s+/).map(Number))
}

function gauss_elimination(A,b){
  const n = A.length
  const M = A.map((row,i)=>row.concat([b[i]]))
  for(let k=0;k<n-1;k++){
    for(let i=k+1;i<n;i++){
      const factor = M[i][k]/M[k][k]
      for(let j=k;j<=n;j++) M[i][j] -= factor*M[k][j]
    }
  }
  const x = Array(n).fill(0)
  for(let i=n-1;i>=0;i--){
    let s = M[i][n]
    for(let j=i+1;j<n;j++) s -= M[i][j]*x[j]
    x[i] = s / M[i][i]
  }
  return x
}

function pivot_partial(A,b){
  const n = A.length
  const M = A.map((row,i)=>row.concat([b[i]]))
  for(let k=0;k<n-1;k++){
    // find max row
    let maxr = k
    for(let i=k+1;i<n;i++){ if(Math.abs(M[i][k])>Math.abs(M[maxr][k])) maxr=i }
    if(maxr!==k){ const tmp = M[k]; M[k]=M[maxr]; M[maxr]=tmp }
    for(let i=k+1;i<n;i++){
      const factor = M[i][k]/M[k][k]
      for(let j=k;j<=n;j++) M[i][j] -= factor*M[k][j]
    }
  }
  const x = Array(n).fill(0)
  for(let i=n-1;i>=0;i--){
    let s=M[i][n]
    for(let j=i+1;j<n;j++) s -= M[i][j]*x[j]
    x[i] = s / M[i][i]
  }
  return x
}

function pivot_total(A,b){
  const n=A.length
  // copy matrix
  const M=A.map((row,i)=>row.concat([b[i]]))
  const marks = Array.from({length:n},(_,i)=>i)
  for(let k=0;k<n-1;k++){
    // find max in submatrix
    let maxv=0, maxr=k, maxc=k
    for(let i=k;i<n;i++) for(let j=k;j<n;j++){
      if(Math.abs(M[i][j])>maxv){ maxv=Math.abs(M[i][j]); maxr=i; maxc=j }
    }
    if(maxv===0) continue
    if(maxr!==k){ const tmp=M[k]; M[k]=M[maxr]; M[maxr]=tmp }
    if(maxc!==k){ for(let i=0;i<n;i++){ const t=M[i][k]; M[i][k]=M[i][maxc]; M[i][maxc]=t } const tmark=marks[k]; marks[k]=marks[maxc]; marks[maxc]=tmark }
    for(let i=k+1;i<n;i++){
      const factor = M[i][k]/M[k][k]
      for(let j=k;j<=n;j++) M[i][j] -= factor*M[k][j]
    }
  }
  const x = Array(n).fill(0)
  for(let i=n-1;i>=0;i--){
    let s=M[i][n]
    for(let j=i+1;j<n;j++) s -= M[i][j]*x[j]
    x[i] = s / M[i][i]
  }
  // reorder solution
  const xfinal = Array(n).fill(0)
  for(let i=0;i<n;i++) xfinal[marks[i]] = x[i]
  return xfinal
}

function renderTable(headers, rows){
  const wrap = $('#tablewrap')
  let html = '<table><thead><tr>' + headers.map(h=>`<th>${h}</th>`).join('') + '</tr></thead><tbody>'
  for(const r of rows){ html += '<tr>' + r.map(c=>`<td>${(typeof c==='number'? c.toPrecision? Number(c).toPrecision(8):c : (c===null? '': c))}</td>`).join('') + '</tr>' }
  html += '</tbody></table>'
  wrap.innerHTML = html
}

$('#run').addEventListener('click', ()=>{
  $('#result').textContent = 'Calculating...'
  $('#tablewrap').innerHTML = ''
  try{
    if(['bisection','false_position','incremental_search','fixed_point','newton','secant','newton_multiple'].includes(current)){
      const expr = $('#fx').value.trim()
      const f = parseF(expr)

      if(current==='incremental_search'){
        const x0 = parseFloat($('#p_x0').value)
        const delta = parseFloat($('#p_delta').value)
        const Nmax = parseInt($('#p_Nmax').value)
        const res = incremental_search(f,x0,delta,Nmax)
        $('#result').textContent = res.interval ? `Interval with sign change: [${res.interval[0]}, ${res.interval[1]}]` : 'No interval found in Nmax iterations'
        renderTable(['i','x0','f(x0)','x1','f(x1)'], res.it)
      }

      else if(current==='bisection'){
        const a=parseFloat($('#p_a').value), b=parseFloat($('#p_b').value)
        const tol=parseFloat($('#p_tol').value), Nmax=parseInt($('#p_Nmax').value)
        const res = bisection(f,a,b,tol,Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','a','b','c','f(c)'], res.it)
      }

      else if(current==='false_position'){
        const a=parseFloat($('#p_a').value), b=parseFloat($('#p_b').value)
        const tol=parseFloat($('#p_tol').value), Nmax=parseInt($('#p_Nmax').value)
        const res = false_position(f,a,b,tol,Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','a','b','c','f(c)'], res.it)
      }

      else if(current==='fixed_point'){
        const gexpr = $('#p_g').value.trim(); const g = parseF(gexpr)
        const x0=parseFloat($('#p_x0').value); const tol=parseFloat($('#p_tol').value); const Nmax=parseInt($('#p_Nmax').value)
        const res = fixed_point(g,x0,tol,Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','x_old','x_new'], res.it)
      }

      else if(current==='newton'){
        const dfexpr = $('#p_df').value.trim(); const df = dfexpr? parseF(dfexpr): null
        const x0=parseFloat($('#p_x0').value); const tol=parseFloat($('#p_tol').value); const Nmax=parseInt($('#p_Nmax').value)
        const res = newton_method(f, df, x0, tol, Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','x','f(x)','f\' (x)','x_next'], res.it)
      }

      else if(current==='secant'){
        const x0=parseFloat($('#p_x0').value), x1=parseFloat($('#p_x1').value)
        const tol=parseFloat($('#p_tol').value); const Nmax=parseInt($('#p_Nmax').value)
        const res = secant(f,x0,x1,tol,Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','x0','x1','x2','f(x1)'], res.it)
      }

      else if(current==='newton_multiple'){
        const dfexpr = $('#p_df').value.trim(); const d2expr = $('#p_d2f').value.trim()
        const df = dfexpr? parseF(dfexpr): null
        const d2 = d2expr? parseF(d2expr): null
        const x0=parseFloat($('#p_x0').value); const tol=parseFloat($('#p_tol').value); const Nmax=parseInt($('#p_Nmax').value)
        const res = newton_multiple(f, df, d2, x0, tol, Nmax)
        $('#result').textContent = `Approximate root: ${res.root}`
        renderTable(['i','x','f(x)','f\'','f\'\'','x_next'], res.it)
      }

    } else if(['gauss','partial_pivot','total_pivot'].includes(current)){
      const n = parseInt($('#p_n').value)
      const Atext = $('#p_A').value
      const bstr = $('#p_b').value
      const A = parseMatrix(Atext)
      const b = bstr.trim().split(/\s+/).map(Number)
      if(A.length!==n) throw new Error('Number of rows in A does not match n')
      let sol
      if(current==='gauss') sol = gauss_elimination(A,b)
      if(current==='partial_pivot') sol = pivot_partial(A,b)
      if(current==='total_pivot') sol = pivot_total(A,b)
      $('#result').textContent = 'Solution: ' + sol.map(v=>v.toPrecision(8)).join(', ')
      $('#tablewrap').innerHTML = ''
    }
  }catch(e){
    $('#result').textContent = 'Error: ' + e.message
    $('#tablewrap').innerHTML = ''
  }
})