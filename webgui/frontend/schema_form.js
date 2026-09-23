// ── schema_form.js — auto-generated run-configuration form ───────────────────
// Fetches /api/schema (field descriptors) and /api/config (current RunSpec as a
// TOML-shaped object), renders a grouped form with NO hard-coded fields, and
// POSTs the edited config back.  Array-of-tables groups (meteorology, species,
// release) are repeatable: each element is its own card with ＋add / ✕remove.

let SCHEMA = [];
let CONFIG = {};

// ── dotted get/set (no array handling) ───────────────────────────────────────
function deepGet(o, dotted) {
  if (!dotted) return o;
  let n = o;
  for (const k of dotted.split('.')) { if (n == null) return undefined; n = n[k]; }
  return n;
}
function deepSet(o, dotted, val) {
  const ks = dotted.split('.');
  let n = o;
  for (let i = 0; i < ks.length - 1; i++) {
    if (n[ks[i]] == null || typeof n[ks[i]] !== 'object') n[ks[i]] = {};
    n = n[ks[i]];
  }
  n[ks[ks.length - 1]] = val;
}

// A field path is either scalar ("output.path") or rooted in one array-of-tables
// ("species[].size.dg" → array "species", sub-path "size.dg").  The "[]" is
// always the first segment.
const ARR_RE = /^([^.[]+)\[\]\.?(.*)$/;
function arrayName(path) { const m = path.match(ARR_RE); return m ? m[1] : null; }

function pathGet(obj, path, idx = 0) {
  const m = path.match(ARR_RE);
  if (m) {
    const a = obj[m[1]];
    if (!Array.isArray(a) || a[idx] == null) return undefined;
    return m[2] ? deepGet(a[idx], m[2]) : a[idx];
  }
  return deepGet(obj, path);
}
function pathSet(obj, path, val, idx = 0) {
  const m = path.match(ARR_RE);
  if (m) {
    if (!Array.isArray(obj[m[1]])) obj[m[1]] = [];
    while (obj[m[1]].length <= idx) obj[m[1]].push({});
    if (m[2]) deepSet(obj[m[1]][idx], m[2], val); else obj[m[1]][idx] = val;
  } else {
    deepSet(obj, path, val);
  }
}

// How many instances of an array group to render (≥1).
function arrCount(arr) {
  return Math.max(1, (Array.isArray(CONFIG[arr]) ? CONFIG[arr].length : 0) || 1);
}

// visible_if "a.b==value" / "a.b!=value", evaluated for element `idx`.
function isVisible(field, idx = 0) {
  const cond = field.visible_if;
  if (!cond) return true;
  let m = cond.match(/^(.*)!=(.*)$/);                 // check != before ==
  if (m) {
    let cur = pathGet(CONFIG, m[1].trim(), idx);
    if (cur === undefined || cur === null) cur = '';
    return String(cur) !== m[2].trim();
  }
  m = cond.match(/^(.*)==(.*)$/);
  if (m) return String(pathGet(CONFIG, m[1].trim(), idx)) === m[2].trim();
  return true;
}

function fieldId(path, idx) { return 'f_' + path.replace(/[^\w]/g, '_') + '__' + idx; }

function inputFor(field, idx) {
  let cur = pathGet(CONFIG, field.path, idx);
  if (cur === undefined) cur = field.default;
  let el;
  if (field.type === 'enum') {
    el = document.createElement('select');
    for (const o of field.options) {
      const opt = document.createElement('option');
      opt.value = o; opt.textContent = o; if (String(cur) === o) opt.selected = true;
      el.appendChild(opt);
    }
  } else if (field.type === 'bool') {
    el = document.createElement('input'); el.type = 'checkbox'; el.checked = !!cur;
  } else if (field.type === 'list') {
    el = document.createElement('input'); el.type = 'text';
    el.value = Array.isArray(cur) ? cur.join(', ') : (cur ?? '');
    el.placeholder = 'comma-separated';
  } else {
    el = document.createElement('input'); el.type = 'text'; el.value = cur ?? '';
  }
  el.id = fieldId(field.path, idx);
  el.dataset.path = field.path; el.dataset.idx = idx; el.dataset.ftype = field.type;
  el.addEventListener('change', () => { collectInto(CONFIG); renderForm(); });
  return el;
}

// Fields that hold a filesystem path get a 📁 browse button (folder vs file).
function pathKind(path) {
  if (/from_file$|lut_path$|variables$/.test(path)) return 'file';
  if (/(^|\.)path$|_dir$|\.dir$/.test(path))        return 'dir';
  return null;
}

function fieldRow(f, idx) {
  const row = document.createElement('div'); row.className = 'schema-row';
  const lab = document.createElement('label');
  lab.textContent = f.path.split('.').pop().replace('[]', '');
  row.dataset.tip = tipText(f);
  row.dataset.tipPath = f.path;
  row.appendChild(lab);
  const input = inputFor(f, idx);
  row.appendChild(input);
  const kind = pathKind(f.path);
  if (kind && typeof browseInto === 'function') {
    const b = document.createElement('button');
    b.className = 'browse-btn'; b.type = 'button'; b.textContent = '📁';
    b.title = kind === 'file' ? 'Browse for a file' : 'Browse for a folder';
    b.onclick = () => browseInto(input.id, kind);
    row.appendChild(b);
  } else if (f.unit) {
    const u = document.createElement('span'); u.className = 'schema-unit'; u.textContent = f.unit; row.appendChild(u);
  }
  return row;
}

// ── hover tooltips ────────────────────────────────────────────────────────────
// One floating card, delegated from #schema-form: hovering any option row pops
// the field's description, unit, default and allowed values (from /api/schema).
function tipText(f) {
  const parts = [];
  if (f.help) parts.push(f.help);
  const meta = [];
  if (f.unit) meta.push(`unit: ${f.unit}`);
  const d = f.default;
  const blank = d === '' || d === null || d === undefined ||
                (Array.isArray(d) && !d.length);
  if (!blank) meta.push(`default: ${Array.isArray(d) ? d.join(', ') : d}`);
  if (f.type === 'enum' && f.options && f.options.length)
    meta.push(`one of: ${f.options.join(' | ')}`);
  if (meta.length) parts.push(meta.join('  ·  '));
  return parts.join('\n');
}

let tipEl = null, tipTimer = null;
function showTip(row) {
  if (!tipEl) { tipEl = document.createElement('div'); tipEl.id = 'schema-tip'; document.body.appendChild(tipEl); }
  tipEl.innerHTML = '';
  const p = document.createElement('div'); p.className = 'tip-path';
  p.textContent = row.dataset.tipPath;
  const b = document.createElement('div'); b.className = 'tip-body';
  b.textContent = row.dataset.tip;
  tipEl.appendChild(p); tipEl.appendChild(b);
  tipEl.style.display = 'block';
  tipEl.style.left = '0px'; tipEl.style.top = '0px';        // measure first
  const w = tipEl.offsetWidth, h = tipEl.offsetHeight;
  const r = row.getBoundingClientRect();
  let x = r.left + 8, y = r.bottom + 6;
  if (x + w > window.innerWidth - 8)  x = window.innerWidth - w - 8;
  if (y + h > window.innerHeight - 8) y = r.top - h - 6;    // flip above
  tipEl.style.left = Math.max(4, x) + 'px';
  tipEl.style.top  = Math.max(4, y) + 'px';
}
function hideTip() {
  clearTimeout(tipTimer); tipTimer = null;
  if (tipEl) tipEl.style.display = 'none';
}
function initTips() {
  const root = document.getElementById('schema-form');
  if (!root || root.dataset.tipsInit) return;
  root.dataset.tipsInit = '1';
  root.addEventListener('mouseover', e => {
    const row = e.target.closest('.schema-row');
    if (!row || !row.dataset.tip) { hideTip(); return; }
    clearTimeout(tipTimer);
    tipTimer = setTimeout(() => showTip(row), 200);
  });
  root.addEventListener('mouseout', e => {
    const row = e.target.closest('.schema-row');
    if (row && row.contains(e.relatedTarget)) return;       // still inside the row
    hideTip();
  });
  document.addEventListener('scroll', hideTip, true);
}

function renderForm() {
  const root = document.getElementById('schema-form');
  if (!root) return;
  root.innerHTML = '';
  const order = []; const byGroup = {};
  for (const f of SCHEMA) { if (!byGroup[f.group]) { byGroup[f.group] = []; order.push(f.group); } byGroup[f.group].push(f); }
  for (const group of order) {
    const fields = byGroup[group];
    const arr = fields.map(f => arrayName(f.path)).find(a => a);   // array name if this is an array group
    const sec = document.createElement('fieldset'); sec.className = 'schema-group';
    const lg = document.createElement('legend'); lg.textContent = group; sec.appendChild(lg);
    if (arr) {
      const count = arrCount(arr);
      for (let i = 0; i < count; i++) {
        const card = document.createElement('div'); card.className = 'arr-elem';
        const hdr = document.createElement('div'); hdr.className = 'arr-elem-hdr';
        const tag = document.createElement('span'); tag.textContent = `#${i + 1}`; hdr.appendChild(tag);
        if (count > 1) {
          const rm = document.createElement('button');
          rm.className = 'browse-btn'; rm.type = 'button'; rm.textContent = '✕'; rm.title = 'Remove';
          rm.onclick = () => removeArrayElem(arr, i);
          hdr.appendChild(rm);
        }
        card.appendChild(hdr);
        for (const f of fields) { if (isVisible(f, i)) card.appendChild(fieldRow(f, i)); }
        sec.appendChild(card);
      }
      const add = document.createElement('button');
      add.className = 'add-btn'; add.type = 'button'; add.textContent = '＋ add ' + arr;
      add.onclick = () => addArrayElem(arr);
      sec.appendChild(add);
    } else {
      for (const f of fields) { if (isVisible(f, 0)) sec.appendChild(fieldRow(f, 0)); }
    }
    root.appendChild(sec);
  }
}

// Build a fresh array element pre-filled with the schema defaults for that group
// (so a new card shows sensible values), skipping blank/empty defaults.
function defaultElem(arr) {
  const e = {};
  for (const f of SCHEMA) {
    const m = f.path.match(ARR_RE);
    if (!m || m[1] !== arr || !m[2]) continue;
    const d = f.default;
    const blank = d === '' || d === undefined || d === null ||
                  (Array.isArray(d) && d.every(x => x === '' || x === undefined || x === null));
    if (!blank) deepSet(e, m[2], d);
  }
  return e;
}
function addArrayElem(arr) {
  collectInto(CONFIG);
  if (!Array.isArray(CONFIG[arr])) CONFIG[arr] = [];
  CONFIG[arr].push(defaultElem(arr));
  renderForm();
}
function removeArrayElem(arr, idx) {
  collectInto(CONFIG);
  if (Array.isArray(CONFIG[arr])) CONFIG[arr].splice(idx, 1);
  if (!CONFIG[arr] || !CONFIG[arr].length) delete CONFIG[arr];
  renderForm();
}

function collectInto(cfg) {
  for (const f of SCHEMA) {
    const arr = arrayName(f.path);
    const count = arr ? arrCount(arr) : 1;
    for (let i = 0; i < count; i++) {
      const el = document.getElementById(fieldId(f.path, i));
      if (!el) continue;
      let v;
      if (f.type === 'bool') v = el.checked;
      else if (f.type === 'int') v = parseInt(el.value, 10);
      else if (f.type === 'number') v = parseFloat(el.value);
      else if (f.type === 'list') {
        v = el.value.trim() === '' ? [] : el.value.split(',').map(s => {
          const t = s.trim(); const n = Number(t); return (t !== '' && !isNaN(n)) ? n : t;
        });
      } else {
        v = el.value;
        if (typeof v === 'string' && v.trim() === '') continue;   // blank = leave unset
      }
      if (v !== undefined && !(typeof v === 'number' && isNaN(v))) pathSet(cfg, f.path, v, i);
    }
  }
}

async function loadSchemaForm() {
  SCHEMA = await (await fetch('/api/schema')).json();
  try { CONFIG = await (await fetch('/api/config')).json(); } catch { CONFIG = {}; }
  if (!CONFIG || Object.keys(CONFIG).length === 0) CONFIG = {};
  renderForm();
  initTips();
}

async function saveSchemaForm() {
  collectInto(CONFIG);
  const res = await fetch('/api/config', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(CONFIG)
  });
  const j = await res.json();
  alert(j.ok ? 'Configuration saved.' : 'Error: ' + (j.msg || 'unknown'));
  return j.ok;
}

// Load a TOML config file from disk (server-side path) into the form.
async function loadTomlConfig() {
  const inp = document.getElementById('toml-path');
  const path = (inp ? inp.value : '').trim();
  if (!path) { alert('Enter the path to a .toml config file.'); return; }
  const res = await fetch('/api/load_toml', {
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ path })
  });
  const j = await res.json();
  if (j.ok === false) { alert('Load error: ' + (j.msg || 'unknown')); return; }
  CONFIG = j;                 // server returns the parsed RunSpec dict
  renderForm();
  if (typeof appendLog === 'function') appendLog('✓ Loaded config: ' + path);
}

// Download the current (edited) configuration as a .toml file.
async function downloadTomlConfig() {
  collectInto(CONFIG);
  await fetch('/api/config', {            // persist first so the server serialises the latest
    method: 'POST', headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(CONFIG)
  });
  const res = await fetch('/api/config.toml');
  if (!res.ok) { alert('Download failed: ' + (await res.text())); return; }
  const text = await res.text();
  const blob = new Blob([text], { type: 'application/toml' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = (pathGet(CONFIG, 'name') || 'easlee_run') + '.toml';
  document.body.appendChild(a); a.click(); a.remove();
  URL.revokeObjectURL(a.href);
}

window.loadSchemaForm     = loadSchemaForm;
window.saveSchemaForm     = saveSchemaForm;
window.loadTomlConfig     = loadTomlConfig;
window.downloadTomlConfig = downloadTomlConfig;
document.addEventListener('DOMContentLoaded', loadSchemaForm);
