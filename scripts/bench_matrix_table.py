# Paper table from the benchmark matrix (scripts/bench_matrix.sh).
# Reads the per-run logs copied into the matrix output directory and prints two markdown
# tables: (A) every row with its cost-centre split, (B) per-process time and speedup at
# 1 / 16 / 64 threads and on the GPU, plus a log-log scaling exponent per process fitted
# on the no-prefetch CPU rows. Groups (timeloop.jl): metread = GRIB read + field
# processing; advection = advance_dispatch! (advection + turbulence, on the GPU incl. the
# host<->device copies); output = deposit_to_grid! (the output KERNEL only); release =
# release + initialize_new_parcels!. "untimed" = TOTAL(run) - sum(groups): the NetCDF
# frame dump, nest reassignment, sync and loop overhead. Convection and deposition are OFF
# in the ETEX protocol and do not appear.
import sys, re, glob, os, math
D = sys.argv[1] if len(sys.argv) > 1 else sorted(glob.glob('/data/jbrioude/easlee_results/bench_matrix_*'))[-1]
G = ['metread', 'advection', 'output', 'release']
rows = {}
for f in sorted(glob.glob(f'{D}/*.log')):
    tag = os.path.basename(f)[:-4]; s = open(f, errors='replace').read()
    m = re.search(r'mode = :([a-z_]+)', s); tot = re.search(r'^\s*TOTAL\(run\)\s+([\d.]+)', s, re.M)
    if not (m and tot): continue
    r = {'mode': m.group(1), 'total': float(tot.group(1))}
    for g in G:
        mm = re.search(rf'^\s*{g}\s+([\d.]+)', s, re.M); r[g] = float(mm.group(1)) if mm else float('nan')
    r['untimed'] = r['total'] - sum(r[g] for g in G)
    mt = re.match(r'(cpu|gpu)_t(\d+)(i?)(_pf)?$', tag)
    r['thr'] = int(mt.group(2)); r['pf'] = bool(mt.group(4)); r['gpu'] = mt.group(1) == 'gpu'
    rows[tag] = r
order = ['cpu_t1','cpu_t2','cpu_t4','cpu_t8','cpu_t16','cpu_t32','cpu_t64','cpu_t16_pf','cpu_t32_pf','cpu_t64_pf','gpu_t4','gpu_t16i','gpu_t4_pf','gpu_t16i_pf']
have = [t for t in order if t in rows]
ref1 = rows.get('cpu_t1'); ref16 = rows.get('cpu_t16')
def f(x, w=7): return f'{x:{w}.1f}' if isinstance(x, float) and math.isfinite(x) else ' ' * (w - 1) + '–'
def sp(a, b): return f'{b/a:.2f}×' if (a and b and math.isfinite(a) and math.isfinite(b) and a > 0) else '–'
print(f'### Table A — every row (run-only seconds; ETEX 1 M parcels, 2 days; `{os.path.basename(D)}`)\n')
print('| row | backend | threads | prefetch | TOTAL | metread | advection | output kernel | release | untimed | vs 1 t | vs 16 t |')
print('|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|')
for t in have:
    r = rows[t]
    print(f"| {t} | {r['mode']} | {r['thr']} | {'on' if r['pf'] else 'off'} | {r['total']:.0f} | {r['metread']:.1f} | {r['advection']:.1f} | {r['output']:.1f} | {r['release']:.1f} | {r['untimed']:.1f} | {sp(r['total'], ref1['total']) if ref1 else '–'} | {sp(r['total'], ref16['total']) if ref16 else '–'} |")
print('\n### Table B — per process: time [s] and speedup (vs 1 thread; in brackets vs 16 threads)\n')
cols = [('cpu_t1', '1 t'), ('cpu_t16', '16 t'), ('cpu_t64', '64 t'), ('gpu_t4', 'GPU (4 t)'), ('gpu_t4_pf', 'GPU + prefetch')]
cols = [c for c in cols if c[0] in rows]
print('| process | ' + ' | '.join(c[1] for c in cols) + ' | scaling exponent (CPU, no prefetch) |')
print('|---|' + '---:|' * (len(cols) + 1))
cpu = [(rows[t]['thr'], rows[t]) for t in have if not rows[t]['gpu'] and not rows[t]['pf']]
def expo(g):
    pts = [(math.log(n), math.log(r[g])) for n, r in cpu if math.isfinite(r[g]) and r[g] > 0]
    if len(pts) < 3: return '–'
    mx = sum(p[0] for p in pts) / len(pts); my = sum(p[1] for p in pts) / len(pts)
    b = sum((p[0]-mx)*(p[1]-my) for p in pts) / sum((p[0]-mx)**2 for p in pts)
    return f't ∝ N^{b:+.2f} ({len(pts)} pts, {min(n for n,_ in cpu)}–{max(n for n,_ in cpu)} t)'
labels = {'metread': 'met fields read + process', 'advection': 'advection + turbulence', 'output': 'output kernel (deposit_to_grid!)', 'release': 'release', 'untimed': 'untimed (frame dump, reassign, sync)', 'total': 'TOTAL(run)'}
for g in ['metread', 'advection', 'output', 'release', 'untimed', 'total']:
    cells = []
    for t, _ in cols:
        v = rows[t][g]; s = f'{v:.1f}'
        if ref1 and t != 'cpu_t1': s += f' ({sp(v, ref1[g])}'
        if ref16 and t not in ('cpu_t1', 'cpu_t16'): s += f'; [{sp(v, ref16[g])}]'
        if ref1 and t != 'cpu_t1': s += ')'
        cells.append(s)
    print(f'| {labels[g]} | ' + ' | '.join(cells) + f' | {expo(g)} |')
print('| convection | off in the ETEX protocol — see the full-physics table (2026-07-09): 4.7× at 6 threads |' + ' |' * len(cols))
print('| dry / wet deposition | off in the ETEX protocol — full-physics table: wet dep 1.5× at 6 threads |' + ' |' * len(cols))
