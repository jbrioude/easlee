// ── State ─────────────────────────────────────────────────────────────────────
let currentView  = 'column';
let lastLogLen   = 0;
let replHistory  = [];
let replHistIdx  = -1;

// ── Tabs ──────────────────────────────────────────────────────────────────────
function switchTab(name) {
  document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
  document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
  document.querySelector(`.tab[data-tab="${name}"]`).classList.add('active');
  document.getElementById(`tab-${name}`).classList.add('active');
}

// ── Status poll (1 s) ─────────────────────────────────────────────────────────
async function pollState() {
  try {
    const res  = await fetch('/api/state');
    const s    = await res.json();

    // Status bar
    const dot  = document.getElementById('status-dot');
    dot.style.color = s.running ? '#a6e3a1' : '#585b70';
    document.getElementById('status-text').textContent     = s.running ? 'running' : 'idle';
    document.getElementById('status-progress').textContent = Math.round(s.progress * 100) + '%';
    document.getElementById('status-time').textContent     = s.time;
    document.getElementById('status-nactive').textContent  = 'active: ' + s.nactive;

    // Progress bar
    const pct = Math.round(s.progress * 100);
    document.getElementById('progress-bar').style.width   = pct + '%';
    document.getElementById('progress-label').textContent = pct + '%';

    // Append new log lines
    if (s.log && s.log.length > lastLogLen) {
      const el = document.getElementById('log-output');
      for (let i = lastLogLen; i < s.log.length; i++) {
        const line = document.createElement('div');
        line.className   = 'log-line';
        line.textContent = s.log[i];
        el.appendChild(line);
      }
      el.scrollTop = el.scrollHeight;
      lastLogLen = s.log.length;
    }

    // Auto-refresh figure when model is running
    if (s.running && document.getElementById('auto-refresh').checked) {
      updateFigure();
    }
    // On run completion, refresh the species list (new output file) + the figure.
    if (wasRunning && !s.running) { refreshSpecies(); updateFigure(); }
    wasRunning = s.running;
  } catch (_) { /* server not ready */ }
}
let wasRunning = false;

// ── Figures ───────────────────────────────────────────────────────────────────
function setView(v) {
  currentView = v;
  if (v === 'column' || v === 'cross') refreshSpecies();
  updateFigure();
}

function updateFigure() {
  const field    = document.getElementById('met-field-sel').value;
  const level    = document.getElementById('level-input').value;
  const lat_frac = document.getElementById('lat-slice').value;
  const species  = document.getElementById('species-sel').value;
  const url = `/api/figure?view=${currentView}&met_field=${field}`
            + `&level=${level}&lat_frac=${lat_frac}`
            + `&species=${encodeURIComponent(species)}&_t=${Date.now()}`;
  document.getElementById('fig-img').src = url;
}

// Populate the species dropdown from the latest output file, keeping the choice.
async function refreshSpecies() {
  try {
    const r = await fetch('/api/species');
    const d = await r.json();
    const sel = document.getElementById('species-sel');
    const cur = sel.value;
    sel.innerHTML = '<option value="">(default)</option>'
      + (d.species || []).map(s => `<option value="${s}">${s}</option>`).join('');
    if ([...sel.options].some(o => o.value === cur)) sel.value = cur;
  } catch (_) { /* server not ready */ }
}

// ── Model control ─────────────────────────────────────────────────────────────
async function controlAction(action) {
  const res  = await fetch('/api/control', {
    method:  'POST',
    headers: {'Content-Type': 'application/json'},
    body:    JSON.stringify({action}),
  });
  const data = await res.json();
  if (!data.ok) appendLog('✗ ' + (data.msg || 'unknown error'));
}

// ── REPL ──────────────────────────────────────────────────────────────────────
function appendReplLine(text, cls) {
  const out  = document.getElementById('repl-output');
  const line = document.createElement('div');
  line.className   = 'repl-line ' + (cls || '');
  line.textContent = text;
  out.appendChild(line);
  document.getElementById('repl-scroll').scrollTop = 99999;
}

async function evalRepl() {
  const input = document.getElementById('repl-input');
  const code  = input.value;
  if (!code.trim()) return;

  appendReplLine('julia> ' + code, 'repl-cmd');
  replHistory.unshift(code);
  replHistIdx = -1;
  input.value = '';

  try {
    const res  = await fetch('/api/eval', {
      method:  'POST',
      headers: {'Content-Type': 'application/json'},
      body:    JSON.stringify({code}),
    });
    const data = await res.json();
    if (data.output) {
      appendReplLine(data.output, data.error ? 'repl-err' : 'repl-ok');
    }
  } catch (e) {
    appendReplLine('Request error: ' + e, 'repl-err');
  }
}

// ── Path picker modal ──────────────────────────────────────────────────────────
let pathTarget = null;   // id of the input being filled
let pathCwd    = '';     // directory currently shown

async function browseInto(inputId, mode = 'any') {
  pathTarget = inputId;
  const start = (document.getElementById(inputId)?.value || '').trim();
  document.getElementById('path-modal-title').textContent =
    mode === 'file' ? 'Select a file' : (mode === 'dir' ? 'Select a folder' : 'Select a folder or file');
  document.getElementById('path-modal').style.display = 'flex';
  await loadPathDir(start);
}
function closePathModal() {
  document.getElementById('path-modal').style.display = 'none';
  pathTarget = null;
}
async function loadPathDir(path) {
  let j;
  try { j = await (await fetch('/api/browse?path=' + encodeURIComponent(path || ''))).json(); }
  catch (_) { return; }
  if (j.ok === false) { if (path) loadPathDir(''); return; }   // bad path → fall back to home
  pathCwd = j.path;
  const sep = j.sep || '/';
  document.getElementById('path-breadcrumb').textContent = j.path;
  document.getElementById('path-current').value = j.path;
  const list = document.getElementById('path-list');
  list.innerHTML = '';
  if (j.parent) list.appendChild(mkPathItem('⬆', '..', true, () => loadPathDir(j.parent)));
  for (const e of j.entries) {
    const full = (j.path.endsWith(sep) ? j.path : j.path + sep) + e.name;
    list.appendChild(mkPathItem(e.dir ? '📁' : '📄', e.name, e.dir,
      e.dir ? () => loadPathDir(full) : () => pickPath(full)));
  }
}
function mkPathItem(ico, name, isDir, onClick) {
  const d = document.createElement('div');
  d.className = 'path-item' + (isDir ? ' is-dir' : '');
  d.innerHTML = `<span class="ico">${ico}</span><span>${escHtml(name)}</span>`;
  d.onclick = onClick;
  return d;
}
function pickPath(p)         { setPathValue(p); closePathModal(); }
function choosePathFolder()  { setPathValue(pathCwd); closePathModal(); }
function setPathValue(p) {
  if (!pathTarget) return;
  const el = document.getElementById(pathTarget);
  if (!el) return;
  el.value = p;
  el.dispatchEvent(new Event('change'));   // let the schema form sync this into the config
}

// ── Release-location map ────────────────────────────────────────────────────────
let mapState = null;   // {lon0,lon1,lat0,lat1,W,H,clicked}
let COAST = null;      // cached coastline polylines

async function ensureCoast() {
  if (COAST !== null) return;
  try { COAST = (await (await fetch('/coastline.json')).json()).lines || []; }
  catch (_) { COAST = []; }
}

async function openMapModal() {
  document.getElementById('map-modal').style.display = 'flex';
  document.getElementById('map-readout').textContent = '';
  document.getElementById('map-status').textContent  = '';
  document.getElementById('map-body').innerHTML =
    '<div style="padding:30px;color:#888">Reading met domain… (first time builds the runtime — may take a moment)</div>';
  await ensureCoast();
  let j;
  try { j = await (await fetch('/api/domain')).json(); }
  catch (_) { document.getElementById('map-body').innerHTML = '<div style="padding:20px;color:#f38ba8">domain request failed</div>'; return; }
  if (j.ok === false || !j.nests || !j.nests.length) {
    document.getElementById('map-body').innerHTML = '<div style="padding:20px;color:#f38ba8">' + (j.msg || 'no met domain available') + '</div>';
    return;
  }
  renderMap(j.nests);
}
function closeMapModal() { document.getElementById('map-modal').style.display = 'none'; }

function niceStep(span) {
  const raw = span / 6, p = Math.pow(10, Math.floor(Math.log10(raw))), m = raw / p;
  const s = m < 1.5 ? 1 : m < 3 ? 2 : m < 7 ? 5 : 10;
  return Math.max(s * p, 0.001);
}
function renderMap(nests) {
  let lon0 = Infinity, lon1 = -Infinity, lat0 = Infinity, lat1 = -Infinity;
  for (const n of nests) {
    lon0 = Math.min(lon0, n.lon_range[0]); lon1 = Math.max(lon1, n.lon_range[1]);
    lat0 = Math.min(lat0, n.lat_range[0]); lat1 = Math.max(lat1, n.lat_range[1]);
  }
  const dlon = (lon1 - lon0) || 1, dlat = (lat1 - lat0) || 1;
  lon0 -= dlon * 0.06; lon1 += dlon * 0.06; lat0 -= dlat * 0.06; lat1 += dlat * 0.06;
  const W = 780;
  const aspect = (lat1 - lat0) / ((lon1 - lon0) * Math.cos((lat0 + lat1) / 2 * Math.PI / 180));
  const H = Math.max(220, Math.min(540, Math.round(W * aspect)));
  mapState = { lon0, lon1, lat0, lat1, W, H, clicked: null };
  const X = lon => (lon - lon0) / (lon1 - lon0) * W;
  const Y = lat => (lat1 - lat) / (lat1 - lat0) * H;
  let s = `<svg id="map-svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" style="background:#0c1320;border-radius:6px;cursor:crosshair;max-width:100%">`;
  const sLon = niceStep(lon1 - lon0), sLat = niceStep(lat1 - lat0);
  for (let lo = Math.ceil(lon0 / sLon) * sLon; lo <= lon1; lo += sLon) {
    const x = X(lo);
    s += `<line x1="${x}" y1="0" x2="${x}" y2="${H}" stroke="#22304a" stroke-width="1"/>`
       + `<text x="${x + 2}" y="${H - 3}" fill="#5b6b85" font-size="10">${lo.toFixed(0)}°</text>`;
  }
  for (let la = Math.ceil(lat0 / sLat) * sLat; la <= lat1; la += sLat) {
    const y = Y(la);
    s += `<line x1="0" y1="${y}" x2="${W}" y2="${y}" stroke="#22304a" stroke-width="1"/>`
       + `<text x="2" y="${y - 2}" fill="#5b6b85" font-size="10">${la.toFixed(0)}°</text>`;
  }
  // coastline: prefer the met land-sea-mask contour delivered per nest with the
  // domain (met-grid resolution); fall back to the bundled low-res coastline.json
  let coastLines = [];
  for (const n of nests) if (n.coast && n.coast.length) coastLines = coastLines.concat(n.coast);
  const coastStroke = coastLines.length ? '#55708f' : '#3b4d68';
  if (!coastLines.length) coastLines = COAST;
  for (const line of coastLines) {
    let inView = false;
    for (const p of line) { if (p[0] >= lon0 && p[0] <= lon1 && p[1] >= lat0 && p[1] <= lat1) { inView = true; break; } }
    if (!inView) continue;
    const d = line.map(p => `${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(' ');
    s += `<polyline points="${d}" fill="none" stroke="${coastStroke}" stroke-width="0.9"/>`;
  }
  const cols = ['#89b4fa', '#a6e3a1', '#f9e2af', '#f38ba8'];
  nests.forEach((n, i) => {
    const pts = n.outline.map(p => `${X(p[0]).toFixed(1)},${Y(p[1]).toFixed(1)}`).join(' ');
    s += `<polygon points="${pts}" fill="${cols[i % 4]}22" stroke="${cols[i % 4]}" stroke-width="1.6"/>`;
    s += `<text x="${X(n.lon_range[0]) + 4}" y="${Y(n.lat_range[1]) + 13}" fill="${cols[i % 4]}" font-size="11">${escHtml(n.model)} (nest ${i + 1})</text>`;
  });
  s += `<circle id="map-marker" r="5" fill="#ff5c8a" stroke="#fff" stroke-width="1.5" style="display:none"/></svg>`;
  document.getElementById('map-body').innerHTML = s;
  document.getElementById('map-svg').addEventListener('click', onMapClick);
  document.getElementById('map-status').textContent = 'Click a location, then "Set release here".';
}
function onMapClick(ev) {
  const svg = document.getElementById('map-svg'), r = svg.getBoundingClientRect();
  const px = (ev.clientX - r.left) / r.width * mapState.W;
  const py = (ev.clientY - r.top) / r.height * mapState.H;
  const lon = mapState.lon0 + px / mapState.W * (mapState.lon1 - mapState.lon0);
  const lat = mapState.lat1 - py / mapState.H * (mapState.lat1 - mapState.lat0);
  mapState.clicked = { lon, lat };
  const mk = document.getElementById('map-marker');
  mk.setAttribute('cx', px); mk.setAttribute('cy', py); mk.style.display = '';
  document.getElementById('map-readout').textContent = `lon ${lon.toFixed(3)}   lat ${lat.toFixed(3)}`;
}
function useMapLocation() {
  if (!mapState || !mapState.clicked) { alert('Click a location on the map first.'); return; }
  const { lon, lat } = mapState.clicked;
  setSchemaField('release[].lon', `${lon.toFixed(3)}, ${lon.toFixed(3)}`);
  setSchemaField('release[].lat', `${lat.toFixed(3)}, ${lat.toFixed(3)}`);
  closeMapModal();
}
function setSchemaField(path, value, idx = 0) {
  // schema inputs are id'd f_<path>__<idx>; the map targets the first release
  const el = document.getElementById('f_' + path.replace(/[^\w]/g, '_') + '__' + idx);
  if (el) { el.value = value; el.dispatchEvent(new Event('change')); }
}

// ── Helpers ───────────────────────────────────────────────────────────────────
function appendLog(msg) {
  const el = document.getElementById('log-output');
  if (!el) return;
  const line = document.createElement('div');
  line.className   = 'log-line';
  line.textContent = msg;
  el.appendChild(line);
  el.scrollTop = el.scrollHeight;
}

function escHtml(s) {
  return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;').replace(/"/g,'&quot;');
}

// ── Boot ──────────────────────────────────────────────────────────────────────
document.addEventListener('DOMContentLoaded', () => {
  // REPL keyboard shortcuts
  const inp = document.getElementById('repl-input');
  inp.addEventListener('keydown', e => {
    if (e.key === 'Enter') {
      evalRepl();
    } else if (e.key === 'ArrowUp') {
      replHistIdx = Math.min(replHistIdx + 1, replHistory.length - 1);
      if (replHistory[replHistIdx] !== undefined) inp.value = replHistory[replHistIdx];
      e.preventDefault();
    } else if (e.key === 'ArrowDown') {
      replHistIdx = Math.max(replHistIdx - 1, -1);
      inp.value   = replHistIdx >= 0 ? replHistory[replHistIdx] : '';
      e.preventDefault();
    }
  });

  // Enter in the TOML path field loads the config
  const tp = document.getElementById('toml-path');
  if (tp) tp.addEventListener('keydown', e => {
    if (e.key === 'Enter' && typeof loadTomlConfig === 'function') loadTomlConfig();
  });

  setInterval(pollState, 1000);
  updateFigure();
});
