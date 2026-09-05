/* Zefiro — interfaccia locale.
   Nessuna fisica qui dentro: tutti i numeri arrivano da /api/v1, che a sua volta
   chiama le stesse funzioni della riga di comando. Questo file disegna e basta. */
'use strict';

const API = '/api/v1';
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const fmt = (v, d = 2) => (v === null || v === undefined || !isFinite(v)) ? '—' : v.toFixed(d);
const SVGNS = 'http://www.w3.org/2000/svg';

const state = { schema: null, design: {}, operating: {}, last: null, thermal: null, runs: [], job: null };

/* ---------------------------------------------------------------- rete --- */
async function api(path, opts = {}) {
  const res = await fetch(API + path, {
    headers: { 'Content-Type': 'application/json' }, ...opts,
    body: opts.body ? JSON.stringify(opts.body) : undefined,
  });
  const txt = await res.text();
  let data = null;
  try { data = txt ? JSON.parse(txt) : null; } catch { data = null; }
  if (!res.ok) {
    const e = (data && data.error) || { code: 'ERRORE_RETE', message: res.statusText, details: {} };
    const err = new Error(e.message); err.code = e.code; err.details = e.details || {};
    throw err;
  }
  return data;
}
const payload = () => ({ design: state.design, operating: state.operating });

/* --------------------------------------------------------------- forma --- */
function el(tag, attrs = {}, ...kids) {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') n.className = v;
    else if (k === 'text') n.textContent = v;
    else if (k.startsWith('on')) n.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined) n.setAttribute(k, v);
  }
  kids.flat().forEach(c => c && n.appendChild(typeof c === 'string' ? document.createTextNode(c) : c));
  return n;
}
function svg(tag, attrs = {}, ...kids) {
  const n = document.createElementNS(SVGNS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v !== null && v !== undefined) n.setAttribute(k, v);
  kids.flat().forEach(c => c && n.appendChild(typeof c === 'string' ? document.createTextNode(c) : c));
  return n;
}
const clear = n => { while (n.firstChild) n.removeChild(n.firstChild); return n; };

function field(f, isOperating) {
  const box = el('div', { class: 'fld' });
  const help = el('div', { class: 'hint', text: f.help || '' });
  const qm = el('button', {
    class: 'qm', type: 'button', 'aria-label': 'spiegazione', text: '?',
    onclick: () => box.classList.toggle('open'),
  });
  const num = el('input', {
    type: 'number', step: f.step || 'any', value: f.value ?? '',
    'aria-label': f.label,
  });
  if (isOperating && f.missing) num.classList.add('miss');

  const row = el('div', { class: 'row' });
  let rng = null;
  if (!isOperating) {
    rng = el('input', {
      type: 'range', min: f.min, max: f.max, step: f.step || 'any',
      value: f.value, 'aria-label': f.label, tabindex: '-1',
    });
    row.appendChild(rng);
  }
  row.appendChild(num);

  const commit = (v, from) => {
    let x = parseFloat(v);
    if (!isFinite(x)) { if (isOperating) { delete state.operating[f.name]; } return; }
    if (!isOperating) {
      x = Math.min(f.max, Math.max(f.min, x));
      if (f.integer) x = Math.round(x);
      state.design[f.name] = x;
      if (from !== 'num') num.value = x;
      if (from !== 'rng' && rng) rng.value = x;
    } else {
      state.operating[f.name] = x;
      num.classList.remove('miss');
    }
  };
  num.addEventListener('input', e => commit(e.target.value, 'num'));
  if (rng) rng.addEventListener('input', e => commit(e.target.value, 'rng'));

  box.append(el('label', {},
    el('span', { class: 'nm', text: f.label }),
    el('span', { class: 'un', text: f.unit }), qm), row, help);
  if (f.value !== null && f.value !== undefined) commit(f.value, 'init');
  return box;
}

function buildForm(s) {
  const host = clear($('#form'));
  const groups = {};
  s.design.forEach(f => (groups[f.group] ||= []).push(f));
  const titles = { ciclo: 'Ciclo', camera: 'Camera', iniezione: 'Iniezione', ugello: 'Ugello' };
  for (const [g, list] of Object.entries(groups)) {
    host.appendChild(el('div', { class: 'grp' },
      el('h3', { text: titles[g] || g }), ...list.map(f => field(f, false))));
  }
  host.appendChild(el('div', { class: 'grp' },
    el('h3', { text: 'Punto operativo' }), ...s.operating.map(f => field(f, true))));

  const c = s.constants;
  host.appendChild(el('div', { class: 'grp' },
    el('h3', { text: 'Fissi (da config)' }),
    el('div', { class: 'hint', style: 'display:block',
      text: `ambiente ${fmt(c.p_amb_bar, 2)} bar · aria ${fmt(c.p_air_supply_bar, 1)} bar (minima a fine raffica) · GPL ${fmt(c.p_fuel_supply_bar, 1)} bar · combustibile ${Object.keys(c.fuel).join('+')} · meccanismo ${c.thermo_source}` })));

  const miss = s.operating.filter(f => f.missing).map(f => f.label);
  const b = $('#dataBadge');
  b.className = 'badge ' + (miss.length ? 'warn' : 'good');
  b.lastElementChild.textContent = miss.length
    ? `${miss.length} dati da misurare` : 'dati operativi completi';
  b.title = miss.join(' · ');
  const rb = $('#revBadge');
  rb.lastElementChild.textContent = s.code_rev;
  rb.title = 'revisione del codice inclusa nel run_id';
}

/* --------------------------------------------------------------- avvisi --- */
function alerts(d) {
  const host = clear($('#alerts'));
  if (!d) return;
  if (d.warnings?.length) host.appendChild(el('div', { class: 'msg warn' },
    el('b', { text: 'Il modello sta lavorando ai suoi limiti' }),
    el('ul', {}, ...d.warnings.map(w => el('li', { text: w })))));
  const viol = Object.entries(d.objectives.g).filter(([, v]) => v > 0);
  if (viol.length) host.appendChild(el('div', { class: 'msg err' },
    el('b', { text: 'Vincoli violati — questo progetto non è realizzabile così' }),
    el('ul', {}, ...viol.map(([k, v]) => el('li', { text: `${k}: ${fmt(v, 3)} (deve essere ≤ 0)` })))));
  if (d.manufacturability?.length) host.appendChild(el('div', { class: 'msg info' },
    el('b', { text: 'Fabbricabilità SLM' }),
    el('ul', {}, ...d.manufacturability.map(m => el('li', { text: m })))));
}

function kpis(d) {
  const host = clear($('#kpis'));
  const l = d.l0, dv = d.derived;
  const items = [
    ['Spinta', fmt(l.thrust, 1), 'N'],
    ['Isp totale', fmt(l.Isp_s, 1), 's'],
    ['T fiamma', fmt(l.T_ad, 0), 'K'],
    ['c*', fmt(l.c_star, 0), 'm/s'],
    ['Rapporto d’area', fmt(l.epsilon, 3), ''],
    ['D gola', fmt(dv.D_t_eq * 1e3, 2), 'mm'],
    ['GPL', fmt((l.mdot_fuel_core + l.mdot_fuel_film) * 3600, 1), 'kg/h'],
    ['L*', fmt(dv.L_star, 2), 'm'],
  ];
  items.forEach(([k, v, u]) => host.appendChild(el('div', { class: 'kpi' },
    el('span', { text: k }), el('b', { text: v }, el('u', { text: u ? ' ' + u : '' })))));
  const p = el('div', { class: 'kpi' }, el('span', { text: 'Esito' }),
    el('b', {}, el('span', {
      class: 'pill ' + (d.objectives.feasible ? 'ok' : 'no'),
      text: d.objectives.feasible ? '✓ fattibile' : '✕ non fattibile',
    })));
  host.appendChild(p);
}

/* -------------------------------------------------------------- disegno --- */
function drawSection(sec) {
  const P = sec.profile_mm, st = sec.stations_mm;
  if (!P.length) return;
  const xs = P.map(p => p[0]), rs = P.map(p => p[1]);
  const x0 = Math.min(...xs), x1 = Math.max(...xs), rmax = Math.max(...rs);
  const M = { l: 62, r: 26, t: 26, b: 46 };
  const W = 900, H = Math.max(260, Math.min(430, (2 * rmax) / (x1 - x0) * (W - M.l - M.r) + M.t + M.b));
  const k = Math.min((W - M.l - M.r) / (x1 - x0), (H - M.t - M.b) / (2 * rmax * 1.12));
  const px = x => M.l + (x - x0) * k;
  const cy = M.t + (H - M.t - M.b) / 2;
  const py = r => cy - r * k;

  const s = clear($('#secSvg'));
  s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const d = P.map((p, i) => `${i ? 'L' : 'M'}${px(p[0]).toFixed(2)},${py(p[1]).toFixed(2)}`).join('') + 'Z';
  const dm = P.map((p, i) => `${i ? 'L' : 'M'}${px(p[0]).toFixed(2)},${(cy + p[1] * k).toFixed(2)}`).join('') + 'Z';

  s.appendChild(svg('rect', { x: 0, y: 0, width: W, height: H, fill: 'var(--gas)' }));
  s.appendChild(svg('path', { d, fill: 'var(--metal)', stroke: 'var(--metal-edge)', 'stroke-width': 1.1, 'stroke-linejoin': 'round' }));
  s.appendChild(svg('path', { d: dm, fill: 'var(--metal)', stroke: 'var(--metal-edge)', 'stroke-width': 1.1, 'stroke-linejoin': 'round' }));
  s.appendChild(svg('line', { x1: M.l - 18, y1: cy, x2: W - M.r + 8, y2: cy, stroke: 'var(--accent)', 'stroke-width': 1, 'stroke-dasharray': '10 3 2 3', opacity: .85 }));

  const mark = (x, label, up) => {
    const X = px(x);
    s.appendChild(svg('line', { x1: X, y1: M.t - 8, x2: X, y2: H - M.b + 6, stroke: 'var(--accent)', 'stroke-width': 1, 'stroke-dasharray': '3 3', opacity: .5 }));
    s.appendChild(svg('text', { x: X, y: up ? M.t - 12 : H - M.b + 20, 'text-anchor': 'middle', class: 'tick', fill: 'var(--accent)' }, label));
  };
  mark(0, 'testa', true);
  mark(st.L_c, 'fine camera', false);
  mark(st.x_lip, 'labbro / gola', true);
  if (st.x_tip > st.x_lip) mark(Math.min(st.x_tip, x1), 'apice plug', false);

  const quota = (xa, xb, y, txt) => {
    const A = px(xa), B = px(xb);
    s.appendChild(svg('line', { x1: A, y1: y, x2: B, y2: y, stroke: 'var(--muted)', 'stroke-width': 1 }));
    [A, B].forEach(X => s.appendChild(svg('line', { x1: X, y1: y - 4, x2: X, y2: y + 4, stroke: 'var(--muted)', 'stroke-width': 1 })));
    s.appendChild(svg('text', { x: (A + B) / 2, y: y - 6, 'text-anchor': 'middle', class: 'tick' }, txt));
  };
  quota(0, st.L_c, H - M.b + 34, `L camera ${fmt(st.L_c, 1)} mm`);
  const yR = py(st.R_c);
  s.appendChild(svg('line', { x1: M.l - 10, y1: yR, x2: M.l - 10, y2: cy, stroke: 'var(--muted)', 'stroke-width': 1 }));
  s.appendChild(svg('text', { x: M.l - 14, y: (yR + cy) / 2, 'text-anchor': 'end', class: 'tick' }, `R ${fmt(st.R_c, 1)}`));

  $('#secCap').textContent =
    `Camera anulare Ø${fmt(2 * st.R_c, 1)} mm lunga ${fmt(st.L_c, 1)} mm · labbro Ø${fmt(2 * st.R_lip, 2)} mm · corpo centrale Ø${fmt(2 * st.r_centerbody, 2)} mm · parete ${fmt(st.t_wall, 2)} mm · gola equivalente Ø${fmt(st.D_t_eq, 2)} mm.`;
}

function drawFace(sec) {
  const j = sec.injectors, st = sec.stations_mm;
  const R = st.R_c + st.t_wall;
  const S = 340, C = S / 2, k = (C - 16) / R;
  const s = clear($('#faceSvg'));
  s.setAttribute('viewBox', `0 0 ${S} ${S}`);
  s.setAttribute('style', 'max-width:360px;margin:0 auto');
  s.appendChild(svg('circle', { cx: C, cy: C, r: R * k, fill: 'var(--metal)', stroke: 'var(--metal-edge)', 'stroke-width': 1.2 }));
  s.appendChild(svg('circle', { cx: C, cy: C, r: st.R_c * k, fill: 'none', stroke: 'var(--metal-edge)', 'stroke-width': 1, 'stroke-dasharray': '4 3', opacity: .8 }));
  s.appendChild(svg('circle', { cx: C, cy: C, r: st.r_centerbody * k, fill: 'var(--metal)', stroke: 'var(--metal-edge)', 'stroke-width': 1.2 }));
  const N = j.N, span = 2 * Math.PI / N;
  const hole = (rad, ang, dia, fill, title) => {
    const g = svg('circle', {
      cx: C + rad * k * Math.cos(ang), cy: C + rad * k * Math.sin(ang),
      r: Math.max(1.4, dia / 2 * k), fill, stroke: 'var(--metal-edge)', 'stroke-width': .7,
    }, svg('title', {}, title));
    s.appendChild(g);
  };
  for (let i = 0; i < N; i++) {
    hole(j.R_inj_mm, i * span + .25 * span, j.d_ox_mm, 'var(--s1)', `aria Ø${fmt(j.d_ox_mm, 2)} mm`);
    hole(j.R_inj_mm, i * span + .75 * span, j.d_fuel_mm, 'var(--s2)', `GPL Ø${fmt(j.d_fuel_mm, 2)} mm`);
    if (j.d_film_mm > 0.05) hole(j.R_film_mm, i * span + .5 * span, j.d_film_mm, 'var(--warning)', `film Ø${fmt(j.d_film_mm, 2)} mm`);
  }
  const leg = [['var(--s1)', `aria Ø${fmt(j.d_ox_mm, 2)}`], ['var(--s2)', `GPL Ø${fmt(j.d_fuel_mm, 2)}`]];
  if (j.d_film_mm > 0.05) leg.push(['var(--warning)', `film Ø${fmt(j.d_film_mm, 2)}`]);
  $('#faceCap').innerHTML = `${N} elementi, passo ${fmt(360 / N, 1)}°` +
    (j.swirl_deg > 0 ? `, swirl ${fmt(j.swirl_deg, 0)}°` : '') + '. ' +
    leg.map(([c, t]) => `<b style="color:${c}">●</b> ${t} mm`).join(' · ');
}

/* --------------------------------------------------------------- grafici --- */
const tip = $('#tip');
function showTip(x, y, html) {
  tip.innerHTML = html; tip.style.opacity = '1';
  const r = tip.getBoundingClientRect();
  tip.style.left = Math.min(window.innerWidth - r.width - 10, x + 14) + 'px';
  tip.style.top = Math.max(8, y - r.height - 12) + 'px';
}
const hideTip = () => { tip.style.opacity = '0'; };

function niceTicks(lo, hi, n = 5) {
  const span = hi - lo || 1;
  const raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || 10 * mag;
  const out = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
  return out;
}

/* Grafico a linee con crocevia e tooltip. Serie sottili, griglia recessiva,
   etichette dirette quando le serie sono poche. */
function lineChart(node, opt) {
  const { series, xLabel, yLabel, threshold, unit = '', xUnit = '' } = opt;
  const W = 900, H = 320, M = { l: 62, r: 146, t: 16, b: 42 };
  const s = clear(node); s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const xs = series[0].x;
  const xlo = Math.min(...xs), xhi = Math.max(...xs);
  let ylo = Infinity, yhi = -Infinity;
  series.forEach(se => se.y.forEach(v => { if (v < ylo) ylo = v; if (v > yhi) yhi = v; }));
  if (threshold) { ylo = Math.min(ylo, threshold.value); yhi = Math.max(yhi, threshold.value); }
  const pad = (yhi - ylo) * 0.12 || 1; ylo -= pad; yhi += pad;
  const PX = v => M.l + (v - xlo) / (xhi - xlo || 1) * (W - M.l - M.r);
  const PY = v => H - M.b - (v - ylo) / (yhi - ylo || 1) * (H - M.t - M.b);

  niceTicks(ylo, yhi).forEach(t => {
    s.appendChild(svg('line', { class: 'gridline', x1: M.l, y1: PY(t), x2: W - M.r, y2: PY(t) }));
    s.appendChild(svg('text', { class: 'tick', x: M.l - 8, y: PY(t) + 3.5, 'text-anchor': 'end' }, String(t)));
  });
  niceTicks(xlo, xhi).forEach(t => s.appendChild(
    svg('text', { class: 'tick', x: PX(t), y: H - M.b + 16, 'text-anchor': 'middle' }, String(t))));
  s.appendChild(svg('line', { class: 'axisline', x1: M.l, y1: H - M.b, x2: W - M.r, y2: H - M.b }));
  s.appendChild(svg('text', { class: 'axlabel', x: M.l, y: H - 6 }, `${xLabel}${xUnit ? ' [' + xUnit + ']' : ''}`));
  s.appendChild(svg('text', { class: 'axlabel', x: 12, y: M.t + 10 }, `${yLabel}${unit ? ' [' + unit + ']' : ''}`));

  if (threshold) {
    s.appendChild(svg('line', {
      x1: M.l, y1: PY(threshold.value), x2: W - M.r, y2: PY(threshold.value),
      stroke: 'var(--critical)', 'stroke-width': 1.5, 'stroke-dasharray': '6 4',
    }));
    s.appendChild(svg('text', {
      x: W - M.r - 4, y: PY(threshold.value) - 6, 'text-anchor': 'end',
      class: 'dlabel', fill: 'var(--critical)',
    }, threshold.label));
  }

  series.forEach(se => {
    s.appendChild(svg('path', {
      d: se.y.map((v, i) => `${i ? 'L' : 'M'}${PX(se.x[i]).toFixed(2)},${PY(v).toFixed(2)}`).join(''),
      fill: 'none', stroke: se.color, 'stroke-width': 2,
      'stroke-linejoin': 'round', 'stroke-linecap': 'round',
      'stroke-dasharray': se.dash || null,
    }));
    const last = se.y.length - 1;
    s.appendChild(svg('text', {
      x: PX(se.x[last]) + 8, y: PY(se.y[last]) + 4, class: 'dlabel', fill: se.color,
    }, se.name));
  });

  const cross = svg('line', { x1: 0, y1: M.t, x2: 0, y2: H - M.b, stroke: 'var(--faint)', 'stroke-width': 1, opacity: 0 });
  s.appendChild(cross);
  const dots = series.map(se => {
    const c = svg('circle', { r: 4.5, fill: se.color, stroke: 'var(--surface)', 'stroke-width': 2, opacity: 0 });
    s.appendChild(c); return c;
  });
  const hit = svg('rect', { x: M.l, y: M.t, width: W - M.l - M.r, height: H - M.t - M.b, fill: 'transparent' });
  s.appendChild(hit);
  hit.addEventListener('pointermove', ev => {
    const box = s.getBoundingClientRect();
    const vx = (ev.clientX - box.left) / box.width * W;
    const xv = xlo + (vx - M.l) / (W - M.l - M.r) * (xhi - xlo);
    let i = 0, best = Infinity;
    xs.forEach((v, k) => { const d = Math.abs(v - xv); if (d < best) { best = d; i = k; } });
    cross.setAttribute('x1', PX(xs[i])); cross.setAttribute('x2', PX(xs[i])); cross.setAttribute('opacity', 1);
    series.forEach((se, k) => {
      dots[k].setAttribute('cx', PX(se.x[i])); dots[k].setAttribute('cy', PY(se.y[i]));
      dots[k].setAttribute('opacity', 1);
    });
    showTip(ev.clientX, ev.clientY,
      `<div class="t">${xLabel} ${fmt(xs[i], 2)} ${xUnit}</div>` +
      series.map(se => `<div class="r"><span><i style="background:${se.color}${se.dash ? ';opacity:.55' : ''}"></i>${se.name}</span><b>${fmt(se.y[i], 0)} ${unit}</b></div>`).join(''));
  });
  hit.addEventListener('pointerleave', () => {
    cross.setAttribute('opacity', 0); dots.forEach(d => d.setAttribute('opacity', 0)); hideTip();
  });
}

function legend(node, items) {
  const h = clear(node);
  items.forEach(([c, t, dash]) => h.appendChild(el('b', {},
    el('i', { style: dash
      ? `background:none;height:0;border-top:2px dashed ${c}`
      : `background:${c}` }), t)));
}

/* --------------------------------------------------------------- termica --- */
function drawThermal(d) {
  const host = clear($('#thmMsg'));
  if (d.material_is_placeholder) host.appendChild(el('div', { class: 'msg warn' },
    el('b', { text: 'Proprietà del 316L non misurate' }),
    document.createTextNode('Il calcolo usa valori segnaposto da manuale per il 316L laminato, che non è il 316L da SLM. I numeri servono a decidere se vale la pena fare la CFD, non a dimensionare. Il carico termico viene da Bartz: correlazione empirica, ±30 %.')));

  const gola = d.stations.find(s => s.label === 'gola');
  const cam = d.stations.find(s => s.label === 'camera');
  // Codifica composita: il COLORE dice la zona (gola / camera), il TRATTO dice
  // la strategia (pieno = pozzo termico, tratteggiato = raffreddato). Due serie
  // dello stesso colore sarebbero indistinguibili con il solo colore.
  const series = [
    { name: 'gola', x: d.t, y: gola.sink, color: 'var(--s2)' },
    { name: 'camera', x: d.t, y: cam.sink, color: 'var(--s1)' },
  ];
  if (gola.cooled) series.splice(1, 0,
    { name: 'gola, con acqua', x: d.t, y: gola.cooled, color: 'var(--s2)', dash: '7 4' });
  lineChart($('#thmSvg'), {
    series: series.slice(0, 3), xLabel: 'tempo', xUnit: 's',
    yLabel: 'temperatura faccia calda', unit: 'K',
    threshold: { value: 1673, label: 'fusione del 316L ≈ 1673 K' },
  });
  legend($('#thmLeg'), [['var(--s2)', 'gola'], ['var(--s1)', 'camera'],
    ['var(--muted)', 'tratteggio = raffreddata ad acqua', true],
    ['var(--critical)', 'limite del materiale', true]]);
  $('#thmCap').textContent = `Temperatura di fiamma ${fmt(d.T_ad, 0)} K, durata ${fmt(d.burn_time, 1)} s. ${d.correlation}.`;

  const t = clear($('#thmTab'));
  t.appendChild(el('thead', {}, el('tr', {},
    ...['zona', 'h gas', 'T parete adiab.', 'q iniziale', 'Biot', 'penetrazione', 'T a fine run', 'con acqua']
      .map(x => el('th', { text: x })))));
  const tb = el('tbody');
  d.stations.forEach(s => tb.appendChild(el('tr', {},
    el('td', { text: s.label }),
    el('td', { class: 'n', text: fmt(s.h_gas, 0) + ' W/m²K' }),
    el('td', { class: 'n', text: fmt(s.T_aw, 0) + ' K' }),
    el('td', { class: 'n', text: fmt(s.q0_MW_m2, 2) + ' MW/m²' }),
    el('td', { class: 'n', text: fmt(s.biot, 2) }),
    el('td', { class: 'n', text: fmt(s.penetration_mm, 1) + ' mm' }),
    el('td', { class: 'n', text: fmt(s.sink[s.sink.length - 1], 0) + ' K' }),
    el('td', { class: 'n', text: s.cooled ? fmt(s.cooled[s.cooled.length - 1], 0) + ' K' : '—' }))));
  t.appendChild(tb);
}

/* --------------------------------------------------------------- storico --- */
function drawRuns(data) {
  state.runs = data.rows;
  const rows = data.rows.filter(r => isFinite(r.thrust) && isFinite(r.Isp_s));
  const s = clear($('#scSvg'));
  const W = 900, H = 340, M = { l: 62, r: 24, t: 16, b: 44 };
  s.setAttribute('viewBox', `0 0 ${W} ${H}`);
  if (!rows.length) { $('#scSvg').appendChild(svg('text', { x: W / 2, y: H / 2, 'text-anchor': 'middle', class: 'tick' }, 'nessuna valutazione ancora')); return; }
  const xs = rows.map(r => r.Isp_s), ys = rows.map(r => r.thrust);
  const xlo = Math.min(...xs), xhi = Math.max(...xs), ylo = Math.min(...ys), yhi = Math.max(...ys);
  const px = v => M.l + (v - xlo) / (xhi - xlo || 1) * (W - M.l - M.r);
  const py = v => H - M.b - (v - ylo) / (yhi - ylo || 1) * (H - M.t - M.b);
  niceTicks(ylo, yhi).forEach(t => {
    s.appendChild(svg('line', { class: 'gridline', x1: M.l, y1: py(t), x2: W - M.r, y2: py(t) }));
    s.appendChild(svg('text', { class: 'tick', x: M.l - 8, y: py(t) + 3.5, 'text-anchor': 'end' }, String(t)));
  });
  niceTicks(xlo, xhi).forEach(t => s.appendChild(svg('text', { class: 'tick', x: px(t), y: H - M.b + 16, 'text-anchor': 'middle' }, String(t))));
  s.appendChild(svg('line', { class: 'axisline', x1: M.l, y1: H - M.b, x2: W - M.r, y2: H - M.b }));
  s.appendChild(svg('text', { class: 'axlabel', x: M.l, y: H - 8 }, 'impulso specifico [s]'));
  s.appendChild(svg('text', { class: 'axlabel', x: 12, y: M.t + 10 }, 'spinta [N]'));

  rows.forEach(r => {
    const X = px(r.Isp_s), Y = py(r.thrust);
    const ok = r.feasible === 1;
    const mark = ok
      ? svg('circle', { cx: X, cy: Y, r: 5, fill: 'var(--good)', 'fill-opacity': .75, stroke: 'var(--surface)', 'stroke-width': 1.5 })
      : svg('path', { d: `M${X - 4},${Y - 4}L${X + 4},${Y + 4}M${X + 4},${Y - 4}L${X - 4},${Y + 4}`, stroke: 'var(--critical)', 'stroke-width': 2, 'stroke-linecap': 'round' });
    mark.style.cursor = 'pointer';
    mark.addEventListener('pointerenter', ev => showTip(ev.clientX, ev.clientY,
      `<div class="t mono">${r.run_id.slice(0, 10)}</div>
       <div class="r"><span>spinta</span><b>${fmt(r.thrust, 1)} N</b></div>
       <div class="r"><span>Isp</span><b>${fmt(r.Isp_s, 1)} s</b></div>
       <div class="r"><span>p camera</span><b>${fmt(r.p_c / 1e5, 2)} bar</b></div>
       <div class="r"><span>φ</span><b>${fmt(r.phi_core, 3)}</b></div>
       <div class="r"><span>esito</span><b>${ok ? '✓ fattibile' : '✕ non fattibile'}</b></div>`));
    mark.addEventListener('pointerleave', hideTip);
    s.appendChild(mark);
  });
  legend($('#scLeg'), [['var(--good)', `● fattibile (${rows.filter(r => r.feasible === 1).length})`], ['var(--critical)', `✕ non fattibile (${rows.filter(r => r.feasible !== 1).length})`]]);

  const t = clear($('#histTab'));
  t.appendChild(el('thead', {}, el('tr', {}, ...['run', 'spinta', 'Isp', 'T fiamma', 'p camera', 'φ', 'ε', 'esito'].map(x => el('th', { text: x })))));
  const tb = el('tbody');
  rows.slice(0, 40).forEach(r => tb.appendChild(el('tr', {},
    el('td', { class: 'mono', text: r.run_id.slice(0, 10) }),
    el('td', { class: 'n', text: fmt(r.thrust, 1) }),
    el('td', { class: 'n', text: fmt(r.Isp_s, 1) }),
    el('td', { class: 'n', text: fmt(r.T_ad, 0) }),
    el('td', { class: 'n', text: fmt(r.p_c / 1e5, 2) }),
    el('td', { class: 'n', text: fmt(r.phi_core, 3) }),
    el('td', { class: 'n', text: fmt(r.epsilon, 3) }),
    el('td', {}, el('span', { class: 'pill ' + (r.feasible === 1 ? 'ok' : 'no'), text: r.feasible === 1 ? 'fattibile' : 'no' })))));
  t.appendChild(tb);
}

/* -------------------------------------------------------------- impianto --- */
function drawPlant(p) {
  const k = clear($('#plantKpis'));
  [['Compressore, limite assoluto', fmt(p.compressore.mdot_isotermo_g_s, 1), 'g/s'],
   ['Massa nel serbatoio', fmt(p.serbatoio.massa_iniziale_g, 0), 'g'],
   ['Estraibile (adiabatico)', fmt(p.serbatoio.utilizzabile_adiabatico_g, 0), 'g'],
   ['T serbatoio a fine raffica', fmt(p.serbatoio.T_finale_C, 0), '°C'],
   [`Portata per ${fmt(p.per_durata.secondi, 0)} s`, fmt(p.per_durata.mdot_g_s, 1), 'g/s'],
  ].forEach(([a, b, u]) => k.appendChild(el('div', { class: 'kpi' },
    el('span', { text: a }), el('b', { text: b }, el('u', { text: ' ' + u })))));

  const pts = p.curva.filter(c => c.adiabatico_s !== null);
  lineChart($('#plantSvg'), {
    series: [{ name: 'adiabatico', x: pts.map(c => c.mdot_g_s), y: pts.map(c => c.adiabatico_s), color: 'var(--s1)' },
             { name: 'isotermo', x: pts.map(c => c.mdot_g_s), y: pts.map(c => c.isotermo_s), color: 'var(--s2)' }],
    xLabel: 'portata d’aria', xUnit: 'g/s', yLabel: 'durata', unit: 's',
  });
}

/* ------------------------------------------------------------------- 3D --- */
/* Visualizzatore STL scritto a mano in WebGL: nessuna libreria, nessuna CDN.
   Il progetto deve funzionare offline e in modo riproducibile, e 600 kB di
   dipendenza per far girare un solido non li vale. */
const viewer = (() => {
  const cv = $('#gl');
  let gl = null, prog = null, buf = null, nTri = 0, raf = 0;
  let theta = -0.55, phi = 0.42, dist = 1.85, cx = 0, cy = 0, cz = 0, radius = 1;
  let wire = false, cut = false;

  const VS = `attribute vec3 p;attribute vec3 n;uniform mat4 mvp;uniform mat3 nm;
    varying vec3 vn;varying vec3 vp;
    void main(){vn=nm*n;vp=p;gl_Position=mvp*vec4(p,1.0);}`;
  const FS = `precision mediump float;varying vec3 vn;varying vec3 vp;
    uniform vec3 col;uniform float cutz;uniform float useCut;
    void main(){
      if(useCut>0.5 && vp.z>cutz) discard;
      vec3 N=normalize(vn);
      float k=max(dot(N,normalize(vec3(0.45,0.75,0.55))),0.0);
      float f=max(dot(N,normalize(vec3(-0.5,-0.3,0.4))),0.0);
      float rim=pow(1.0-abs(N.z),2.5);
      vec3 c=col*(0.24+0.68*k+0.20*f)+vec3(0.16)*rim;
      gl_FragColor=vec4(c,1.0);}`;

  function shader(type, src) {
    const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
    return s;
  }
  function init() {
    if (gl) return true;
    gl = cv.getContext('webgl', { antialias: true, alpha: false });
    if (!gl) return false;
    prog = gl.createProgram();
    gl.attachShader(prog, shader(gl.VERTEX_SHADER, VS));
    gl.attachShader(prog, shader(gl.FRAGMENT_SHADER, FS));
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    gl.useProgram(prog);
    gl.enable(gl.DEPTH_TEST);
    buf = gl.createBuffer();
    bindEvents();
    return true;
  }
  function parseSTL(ab) {
    const dv = new DataView(ab);
    const n = dv.getUint32(80, true);
    if (84 + n * 50 !== ab.byteLength) throw new Error('STL non binario o corrotto');
    const data = new Float32Array(n * 18);
    let o = 84, w = 0;
    let x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let i = 0; i < n; i++) {
      const nx = dv.getFloat32(o, true), ny = dv.getFloat32(o + 4, true), nz = dv.getFloat32(o + 8, true);
      o += 12;
      for (let v = 0; v < 3; v++) {
        const X = dv.getFloat32(o, true), Y = dv.getFloat32(o + 4, true), Z = dv.getFloat32(o + 8, true);
        o += 12;
        data[w++] = X; data[w++] = Y; data[w++] = Z;
        data[w++] = nx; data[w++] = ny; data[w++] = nz;
        if (X < x0) x0 = X; if (X > x1) x1 = X;
        if (Y < y0) y0 = Y; if (Y > y1) y1 = Y;
        if (Z < z0) z0 = Z; if (Z > z1) z1 = Z;
      }
      o += 2;
    }
    return { data, n, box: [x0, y0, z0, x1, y1, z1] };
  }
  function load(ab) {
    if (!init()) throw new Error('WebGL non disponibile');
    const { data, n, box } = parseSTL(ab);
    nTri = n;
    cx = (box[0] + box[3]) / 2; cy = (box[1] + box[4]) / 2; cz = (box[2] + box[5]) / 2;
    radius = Math.max(box[3] - box[0], box[4] - box[1], box[5] - box[2]) / 2 || 1;
    dist = 1.85;
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW);
    $('#glInfo').textContent = `${n.toLocaleString('it')} triangoli · ingombro ${(box[3] - box[0]).toFixed(1)} × ${(box[4] - box[1]).toFixed(1)} × ${(box[5] - box[2]).toFixed(1)} mm`;
    draw();
  }
  function mul(a, b) {
    const r = new Float32Array(16);
    for (let i = 0; i < 4; i++) for (let j = 0; j < 4; j++) {
      let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + j] * b[i * 4 + k];
      r[i * 4 + j] = s;
    }
    return r;
  }
  function draw() {
    if (!gl || !nTri) return;
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = cv.clientWidth, h = cv.clientHeight;
    if (cv.width !== w * dpr || cv.height !== h * dpr) { cv.width = w * dpr; cv.height = h * dpr; }
    gl.viewport(0, 0, cv.width, cv.height);
    const bg = getComputedStyle(document.body).getPropertyValue('--surface-2').trim();
    const m = /#(\w\w)(\w\w)(\w\w)/.exec(bg);
    gl.clearColor(m ? parseInt(m[1], 16) / 255 : .93, m ? parseInt(m[2], 16) / 255 : .94, m ? parseInt(m[3], 16) / 255 : .95, 1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

    const s = 1 / radius;
    const T = new Float32Array([s, 0, 0, 0, 0, s, 0, 0, 0, 0, s, 0, -cx * s, -cy * s, -cz * s, 1]);
    const ct = Math.cos(theta), st = Math.sin(theta), cp = Math.cos(phi), sp = Math.sin(phi);
    const Ry = new Float32Array([ct, 0, -st, 0, 0, 1, 0, 0, st, 0, ct, 0, 0, 0, 0, 1]);
    const Rx = new Float32Array([1, 0, 0, 0, 0, cp, sp, 0, 0, -sp, cp, 0, 0, 0, 0, 1]);
    const V = new Float32Array([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, -dist, 1]);
    const asp = w / h, f = 1 / Math.tan(0.62), far = 100, near = 0.05;
    const P = new Float32Array([f / asp, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) / (near - far), -1, 0, 0, 2 * far * near / (near - far), 0]);
    const MV = mul(V, mul(Rx, mul(Ry, T)));
    gl.uniformMatrix4fv(gl.getUniformLocation(prog, 'mvp'), false, mul(P, MV));
    gl.uniformMatrix3fv(gl.getUniformLocation(prog, 'nm'), false, new Float32Array([
      MV[0], MV[1], MV[2], MV[4], MV[5], MV[6], MV[8], MV[9], MV[10]]));
    const dark = getComputedStyle(document.body).getPropertyValue('--ground').trim().length && matchMedia('(prefers-color-scheme: dark)').matches;
    gl.uniform3f(gl.getUniformLocation(prog, 'col'), dark ? .42 : .62, dark ? .49 : .68, dark ? .55 : .73);
    gl.uniform1f(gl.getUniformLocation(prog, 'cutz'), cz);
    gl.uniform1f(gl.getUniformLocation(prog, 'useCut'), cut ? 1 : 0);

    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    const pa = gl.getAttribLocation(prog, 'p'), na = gl.getAttribLocation(prog, 'n');
    gl.enableVertexAttribArray(pa); gl.vertexAttribPointer(pa, 3, gl.FLOAT, false, 24, 0);
    gl.enableVertexAttribArray(na); gl.vertexAttribPointer(na, 3, gl.FLOAT, true, 24, 12);
    gl.drawArrays(gl.TRIANGLES, 0, nTri * 3);
    if (wire) {
      gl.uniform3f(gl.getUniformLocation(prog, 'col'), 0.05, 0.05, 0.06);
      for (let i = 0; i < nTri; i++) gl.drawArrays(gl.LINE_LOOP, i * 3, 3);
    }
  }
  const schedule = () => { if (!raf) raf = requestAnimationFrame(() => { raf = 0; draw(); }); };
  function bindEvents() {
    let drag = false, lx = 0, ly = 0;
    cv.addEventListener('pointerdown', e => { drag = true; lx = e.clientX; ly = e.clientY; cv.setPointerCapture(e.pointerId); });
    cv.addEventListener('pointerup', e => { drag = false; try { cv.releasePointerCapture(e.pointerId); } catch {} });
    cv.addEventListener('pointermove', e => {
      if (!drag) return;
      theta += (e.clientX - lx) * 0.008; phi += (e.clientY - ly) * 0.008;
      phi = Math.max(-1.5, Math.min(1.5, phi)); lx = e.clientX; ly = e.clientY; schedule();
    });
    cv.addEventListener('wheel', e => {
      e.preventDefault();
      dist = Math.max(1.15, Math.min(9, dist * (1 + Math.sign(e.deltaY) * 0.11))); schedule();
    }, { passive: false });
    window.addEventListener('resize', schedule);
    matchMedia('(prefers-color-scheme: dark)').addEventListener('change', schedule);
    $('#glWire').addEventListener('change', e => { wire = e.target.checked; schedule(); });
    $('#glCut').addEventListener('change', e => { cut = e.target.checked; schedule(); });
  }
  return { load, redraw: schedule };
})();

/* ----------------------------------------------------------------- lavori --- */
async function runJob(kind, body, label) {
  const bar = $('#progBar'), msg = $('#jobMsg');
  msg.textContent = label + '…'; bar.style.width = '3%';
  const job = await api(`/jobs/${kind}`, { method: 'POST', body });
  state.job = job.id;
  for (;;) {
    await new Promise(r => setTimeout(r, 350));
    const st = await api(`/jobs/${job.id}`);
    bar.style.width = Math.max(3, st.progress * 100) + '%';
    if (st.message) msg.textContent = st.message;
    if (st.status === 'done') { bar.style.width = '100%'; msg.textContent = label + ' completato'; setTimeout(() => { bar.style.width = '0'; }, 900); return st.result; }
    if (st.status === 'failed') { bar.style.width = '0'; throw Object.assign(new Error(st.error?.message || 'fallito'), { code: st.error?.code }); }
    if (st.status === 'cancelled') { bar.style.width = '0'; msg.textContent = 'annullato'; return null; }
  }
}

function showError(err) {
  const host = $('#alerts');
  host.insertBefore(el('div', { class: 'msg err' },
    el('b', { text: err.code || 'Errore' }),
    document.createTextNode(err.message),
    err.details?.fields ? el('ul', {}, ...err.details.fields.map(f => el('li', { text: f }))) : null,
  ), host.firstChild);
}

/* ------------------------------------------------------------------ flusso --- */
async function evaluate() {
  const btn = $('#btnEval'); btn.disabled = true;
  try {
    const d = await api('/evaluate', { method: 'POST', body: payload() });
    state.last = d;
    alerts(d); kpis(d); drawSection(d.section); drawFace(d.section);
    $('#dlStl').href = `${API}/geometry/${d.run_id}/stl`;
    $('#dlStep').href = `${API}/geometry/${d.run_id}/step`;
    try { drawThermal(state.thermal = await api('/thermal', { method: 'POST', body: payload() })); }
    catch (e) { clear($('#thmMsg')).appendChild(el('div', { class: 'msg warn' }, el('b', { text: 'Termica non disponibile' }), document.createTextNode(e.message))); }
    drawRuns(await api('/runs?limit=500'));
  } catch (e) { clear($('#alerts')); showError(e); }
  finally { btn.disabled = false; }
}

async function geometry() {
  const btn = $('#btnGeom'); btn.disabled = true;
  try {
    const g = await runJob('geometry', { ...payload(), sector: false }, 'Costruzione del solido');
    if (!g) return;
    $('#dlStl').href = g.stl_url; $('#dlStep').href = g.step_url;
    selectTab('v3d');
    viewer.load(await (await fetch(g.stl_url)).arrayBuffer());
    $('#glInfo').textContent += ` · massa ${fmt(g.mass_g, 0)} g (densità segnaposto) · ${g.watertight ? 'chiuso ✓' : 'NON chiuso ✕'}`;
  } catch (e) { showError(e); } finally { btn.disabled = false; }
}

async function sweep() {
  const btn = $('#btnSweep'); btn.disabled = true;
  try {
    const r = await runJob('sweep', { n: 200, seed: 0, operating: state.operating }, 'Piano sperimentale');
    if (r) { $('#jobMsg').textContent = `${r.valutati} valutati, ${r.scartati} scartati, ${r.fattibili_totali} fattibili su ${r.run_totali}`; selectTab('hist'); drawRuns(await api('/runs?limit=1000')); }
  } catch (e) { showError(e); } finally { btn.disabled = false; }
}

function selectTab(v) {
  $$('.tab').forEach(t => t.setAttribute('aria-selected', String(t.dataset.v === v)));
  $$('.view').forEach(s => s.classList.toggle('on', s.id === 'v-' + v));
  if (v === 'v3d') viewer.redraw();
}

async function boot() {
  try {
    state.schema = await api('/schema');
    buildForm(state.schema);
    $$('.tab').forEach(t => t.addEventListener('click', () => selectTab(t.dataset.v)));
    $('#btnEval').addEventListener('click', evaluate);
    $('#btnGeom').addEventListener('click', geometry);
    $('#btnSweep').addEventListener('click', sweep);
    drawPlant(await api('/plant'));
    drawRuns(await api('/runs?limit=500'));

    // Se mancano dati operativi, non si lancia una valutazione destinata a
    // fallire: si dice che cosa manca. Un errore all'avvio farebbe sembrare
    // rotto un software che sta invece facendo il suo lavoro.
    const miss = state.schema.operating.filter(f => f.missing);
    if (miss.length) {
      clear($('#alerts')).appendChild(el('div', { class: 'msg warn' },
        el('b', { text: 'Mancano dei dati misurati' }),
        document.createTextNode('Il punto operativo in config/ ha dei valori non ancora misurati. Compilali qui a sinistra (restano in memoria, non toccano i file) oppure misurali e mettili in config/operating_point.yaml.'),
        el('ul', {}, ...miss.map(f => el('li', { text: `${f.label} — ${f.help}` })))));
      clear($('#kpis'));
    } else {
      await evaluate();
    }
  } catch (e) { showError(e); }
}
boot();
