#!/usr/bin/env python3
"""plot_estimator_2x2.py — figures for ESTIMATOR_2X2.md (the dual-estimator
benchmark, ETEX + CAPTEX R1, 2026-07-06/07).

Re-pairs the four estimator cells from the RAW run outputs (the two grids +
the two receptor CSVs per experiment) against the DATEM observations — the
same pairing as scripts/check_etex.jl / check_captex.jl — recomputes the
ATMES-II statistics (must reproduce ETEX_RESULTS.md / CAPTEX_RESULTS.md),
and writes the explanatory figures to results/estimator_2x2/.

  python3 scripts/plot_estimator_2x2.py
"""
import os, re, math
from datetime import datetime, timedelta

import numpy as np
import netCDF4
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle, Circle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JULIA = os.path.dirname(ROOT)
OUT = os.path.join(ROOT, "results", "estimator_2x2")
os.makedirs(OUT, exist_ok=True)

C_UNI, C_AB, C_REC = "#0072B2", "#D55E00", "#009E73"   # uniform / adaptive / receptor
BAND_EDGES = [0, 150, 400, 800, 2000]
BAND_COLS = ["#c02020", "#e08020", "#4090c0", "#606060"]

# ── DATEM parsing (mirror scripts/datem.jl) ───────────────────────────────────
def read_meas(path):
    obs = []
    for ln in open(path):
        if not re.match(r"^\s*\d{4}\s", ln):
            continue
        f = ln.split()
        if len(f) < 9:
            continue
        yr, mo, dy = int(f[0]), int(f[1]), int(f[2])
        sh, sm = int(f[3][:2]), int(f[3][2:4])
        dh, dm = int(f[4][:2]), int(f[4][2:4])
        start = datetime(yr, mo, dy, sh, sm)
        obs.append(dict(start=start, stop=start + timedelta(hours=dh, minutes=dm),
                        lat=float(f[5]), lon=float(f[6]),
                        conc=float(f[7]), stn=int(f[8])))
    return obs

def datem_stats(o, m):
    o = np.asarray(o, float); m = np.asarray(m, float)
    n = len(o)
    ob, mb = o.mean(), m.mean()
    fb = 2 * (mb - ob) / (mb + ob) if (mb + ob) > 0 else 0.0
    nmse = ((m - o) ** 2).mean() / (mb * ob) if mb * ob > 0 else np.nan
    r = float(np.corrcoef(o, m)[0, 1]) if o.std() > 0 and m.std() > 0 else 0.0
    both = np.sum((o > 0) & (m > 0)); either = np.sum((o > 0) | (m > 0))
    fms = 100.0 * both / either if either > 0 else 0.0
    xs = np.unique(np.concatenate([o, m]))
    co = np.searchsorted(np.sort(o), xs, side="right") / n
    cm = np.searchsorted(np.sort(m), xs, side="right") / n
    ksp = 100.0 * np.abs(co - cm).max()
    rank = r ** 2 + (1 - abs(fb) / 2) + fms / 100 + (1 - ksp / 100)
    return dict(n=n, obar=ob, mbar=mb, fb=fb, nmse=nmse, r=r,
                fms=fms, ksp=ksp, rank=rank)

def haversine_km(lat1, lon1, lat2, lon2):
    p = math.pi / 180
    a = (math.sin((lat2 - lat1) * p / 2) ** 2 +
         math.cos(lat1 * p) * math.cos(lat2 * p) * math.sin((lon2 - lon1) * p / 2) ** 2)
    return 2 * 6371.0 * math.asin(math.sqrt(a))

# ── model output readers ──────────────────────────────────────────────────────
class Grid:
    """Time-framed grid + the same bilinear station sampling as check_*.jl."""
    def __init__(self, path):
        ds = netCDF4.Dataset(path)
        self.lon = ds["lon"][:].filled(np.nan)
        self.lat = ds["lat"][:].filled(np.nan)
        self.conc = ds["conc"][:].filled(0.0)            # (time, height, lat, lon)
        t0 = datetime.fromisoformat(re.search(r"since (.+)$", ds["time"].units).group(1))
        self.frame_end = [t0 + timedelta(seconds=float(s)) for s in ds["time"][:]]
        self.avg_s = int(round(float(ds.average_s)))
        self.t0 = t0
        ds.close()
        self.klay = 1 if self.conc.shape[1] >= 2 else 0
        self.dlon = self.lon[1] - self.lon[0]
        self.dlat = self.lat[1] - self.lat[0]
        self.frame_of = {t: i for i, t in enumerate(self.frame_end)}

    def sample_pg(self, lo, la, ti):
        gx = (lo - self.lon[0]) / self.dlon
        gy = (la - self.lat[0]) / self.dlat
        if not (0 <= gx <= len(self.lon) - 1 and 0 <= gy <= len(self.lat) - 1):
            return np.nan
        i0 = min(max(int(gx), 0), len(self.lon) - 2); wx = gx - i0
        j0 = min(max(int(gy), 0), len(self.lat) - 2); wy = gy - j0
        c = self.conc[ti, self.klay]
        v = ((1 - wx) * (1 - wy) * c[j0, i0] + wx * (1 - wy) * c[j0, i0 + 1] +
             (1 - wx) * wy * c[j0 + 1, i0] + wx * wy * c[j0 + 1, i0 + 1])
        return float(v) * 1e15                            # kg/m³ → pg/m³

    def window_pg(self, ob):
        """3-h window = its frame; 6-h window = mean of its two frames."""
        dur = int((ob["stop"] - ob["start"]).total_seconds())
        if dur == self.avg_s and ob["stop"] in self.frame_of:
            return self.sample_pg(ob["lon"], ob["lat"], self.frame_of[ob["stop"]])
        half = ob["stop"] - timedelta(seconds=self.avg_s)
        if dur == 2 * self.avg_s and ob["stop"] in self.frame_of and half in self.frame_of:
            a = self.sample_pg(ob["lon"], ob["lat"], self.frame_of[half])
            b = self.sample_pg(ob["lon"], ob["lat"], self.frame_of[ob["stop"]])
            return 0.5 * (a + b)
        return np.nan

class Receptors:
    """Gridless KDE dump: time_s,station,conc_kgm3,resolution_m,nsamp."""
    def __init__(self, path, t0, avg_s):
        self.val, self.res = {}, {}
        self.avg_s = avg_s
        for ln in open(path):
            if ln.startswith("#") or ln.startswith("time_s"):
                continue
            f = ln.strip().split(",")
            if len(f) < 5:
                continue
            t = t0 + timedelta(seconds=int(f[0])); stn = int(f[1])
            self.val[(stn, t)] = float(f[2]) * 1e15
            try:
                self.res[(stn, t)] = float(f[3])
            except ValueError:
                pass

    def window_pg(self, ob):
        dur = int((ob["stop"] - ob["start"]).total_seconds())
        k = (ob["stn"], ob["stop"])
        if dur == self.avg_s and k in self.val:
            return self.val[k]
        k2 = (ob["stn"], ob["stop"] - timedelta(seconds=self.avg_s))
        if dur == 2 * self.avg_s and k in self.val and k2 in self.val:
            return 0.5 * (self.val[k] + self.val[k2])
        return np.nan

# ── pair one experiment ───────────────────────────────────────────────────────
def pair_experiment(obs, grids, receptors, src):
    """Rows: one per station-window present in ALL four estimators."""
    rows = []
    for ob in obs:
        vals = {}
        ok = True
        for name, g in grids.items():
            v = g.window_pg(ob)
            if not np.isfinite(v):
                ok = False; break
            vals["grid_" + name] = v
        if not ok:
            continue
        for name, rcp in receptors.items():
            v = rcp.window_pg(ob)
            if not np.isfinite(v):
                ok = False; break
            vals["rec_" + name] = v
        if not ok:
            continue
        rows.append(dict(stn=ob["stn"], lat=ob["lat"], lon=ob["lon"],
                         stop=ob["stop"], obs=ob["conc"],
                         dist=haversine_km(src[0], src[1], ob["lat"], ob["lon"]),
                         **vals))
    return rows

def col(rows, k):
    return np.array([r[k] for r in rows])

# ── load ETEX ─────────────────────────────────────────────────────────────────
ETEX_OBS = os.path.join(JULIA, "ETEX")
ETEX_DIR = os.path.join(ROOT, "results", "etex")
ETEX_SRC = (48.058, -2.0083)                             # Monterfil (emit-t1.txt)

obs_e = read_meas(os.path.join(ETEX_OBS, "meas-t1.txt"))
g_uni_e = Grid(os.path.join(ETEX_DIR, "grid_output_etex_uni.nc"))
g_ab_e  = Grid(os.path.join(ETEX_DIR, "grid_output_etex_ab.nc"))
r_uni_e = Receptors(os.path.join(ETEX_DIR, "receptors_output_etex_uni.csv"), g_uni_e.t0, g_uni_e.avg_s)
r_ab_e  = Receptors(os.path.join(ETEX_DIR, "receptors_output_etex_ab.csv"),  g_ab_e.t0,  g_ab_e.avg_s)
etex = pair_experiment(obs_e, {"uni": g_uni_e, "ab": g_ab_e},
                       {"uni": r_uni_e, "ab": r_ab_e}, ETEX_SRC)

# ── load CAPTEX R1 ────────────────────────────────────────────────────────────
CAP_OBS = os.path.join(JULIA, "CAPTEX")
CAP_DIR = os.path.join(ROOT, "results", "captex")
CAP_SRC = (39.80, -84.05)                                # release 1 (emit-t1.txt)

obs_c = read_meas(os.path.join(CAP_OBS, "meas-t1.txt"))
g_uni_c = Grid(os.path.join(CAP_DIR, "grid_output_captex_uni.nc"))
g_ab_c  = Grid(os.path.join(CAP_DIR, "grid_output_captex_ab.nc"))
r_uni_c = Receptors(os.path.join(CAP_DIR, "receptors_output_captex_uni.csv"), g_uni_c.t0, g_uni_c.avg_s)
r_ab_c  = Receptors(os.path.join(CAP_DIR, "receptors_output_captex_ab.csv"),  g_ab_c.t0,  g_ab_c.avg_s)
captex = pair_experiment(obs_c, {"uni": g_uni_c, "ab": g_ab_c},
                         {"uni": r_uni_c, "ab": r_ab_c}, CAP_SRC)

CELLS = [("grid_uni", "grid · uniform", C_UNI),
         ("grid_ab",  "grid · adaptive", C_AB),
         ("rec_uni",  "receptors (uniform run)", C_REC),
         ("rec_ab",   "receptors (adaptive run)", "#7bc8a4")]

def print_table(tag, rows):
    print(f"\n{tag}  ({len(rows)} station-windows)")
    stats = {}
    for key, label, _ in CELLS:
        s = datem_stats(col(rows, "obs"), col(rows, key))
        stats[key] = s
        print(f"  {label:26s} R={s['r']:.3f} FB={s['fb']:+.2f} "
              f"FMS={s['fms']:.0f}% KSP={s['ksp']:.1f} RANK={s['rank']:.2f}")
    return stats

stats_e = print_table("ETEX (must reproduce ETEX_RESULTS.md 2x2)", etex)
stats_c = print_table("CAPTEX R1 (must reproduce CAPTEX_RESULTS.md 2x2)", captex)

# dump the joined pairs for reuse
for name, rows in (("etex_pairs.csv", etex), ("captex_pairs.csv", captex)):
    with open(os.path.join(OUT, name), "w") as io:
        io.write("stn,lat,lon,window_end,dist_km,obs_pg,grid_uni_pg,grid_ab_pg,rec_uni_pg,rec_ab_pg\n")
        for r in rows:
            io.write(f"{r['stn']},{r['lat']:.3f},{r['lon']:.3f},{r['stop'].isoformat()},"
                     f"{r['dist']:.1f},{r['obs']:.2f},{r['grid_uni']:.2f},{r['grid_ab']:.2f},"
                     f"{r['rec_uni']:.4g},{r['rec_ab']:.4g}\n")

# ═══ FIG 1 — what the four estimators ARE (schematic) ═════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.6))
rng = np.random.default_rng(4)
# a tight filament of parcels passing north of a station (the 8021 geometry)
n_p = 90
fil_x = np.sort(rng.uniform(0.4, 5.6, n_p))
fil_y = 3.62 + 0.35 * np.sin(fil_x * 0.9) + rng.normal(0, 0.06, n_p)
stn_xy = (3.1, 1.9)

def draw_grid(ax, cell=1.0, n=6, lw=0.7):
    for i in range(n + 1):
        ax.plot([0, n], [i, i], color="0.75", lw=lw, zorder=1)
        ax.plot([i, i], [0, n], color="0.75", lw=lw, zorder=1)

def cells_hit(ax, xs, ys, spread, alpha_hi):
    """Shade cells within `spread` cells of any parcel (the deposit footprint)."""
    hit = np.zeros((6, 6))
    for x, y in zip(xs, ys):
        for j in range(6):
            for i in range(6):
                cx, cy = i + 0.5, j + 0.5
                d = max(abs(cx - x), abs(cy - y))
                if d <= spread:
                    hit[j, i] += max(0.0, 1 - d / (spread + 1e-9))
    hit /= hit.max()
    for j in range(6):
        for i in range(6):
            if hit[j, i] > 0.02:
                ax.add_patch(Rectangle((i, j), 1, 1, color="#e08020",
                                       alpha=alpha_hi * hit[j, i], zorder=2, lw=0))
    return hit

for ax in axes:
    ax.set_xlim(0, 6); ax.set_ylim(0, 6)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])

# (a) uniform
ax = axes[0]
draw_grid(ax)
cells_hit(ax, fil_x, fil_y, 1.35, 0.75)
ax.plot(fil_x, fil_y, ".", color="k", ms=3, zorder=4)
ax.plot(*stn_xy, "^", color=C_UNI, ms=13, mec="k", zorder=5)
ax.add_patch(Rectangle((stn_xy[0] - 1.5, stn_xy[1] - 1.5), 3, 3, fill=False,
                       ec=C_UNI, lw=1.8, ls="--", zorder=5))
ax.annotate("bilinear deposit: each parcel\nspreads over ±1 cell",
            xy=(2.8, 5.28), ha="center", fontsize=10)
ax.annotate("bilinear sampling:\nstation reads ±1 cell", xy=(stn_xy[0], 0.3),
            ha="center", fontsize=10, color=C_UNI)
ax.set_title("(a) grid — UNIFORM kernel\nmass smeared twice → station sees the plume",
             fontsize=11)

# (b) adaptive
ax = axes[1]
draw_grid(ax)
cells_hit(ax, fil_x, fil_y, 0.28, 0.9)
ax.plot(fil_x, fil_y, ".", color="k", ms=3, zorder=4)
ax.plot(*stn_xy, "^", color=C_AB, ms=13, mec="k", zorder=5)
for x, y in zip(fil_x[::9], fil_y[::9]):
    ax.add_patch(Circle((x, y), 0.16, fill=False, ec="#c02020", lw=1.0, zorder=3))
ax.annotate("kernel width h = local parcel\nspacing (kNN) ≪ cell size", xy=(2.8, 5.05),
            ha="center", fontsize=10, color="#c02020")
ax.annotate("station cell stays EMPTY\n→ the miss is reported", xy=(stn_xy[0], 0.3),
            ha="center", fontsize=10, color=C_AB)
ax.set_title("(b) grid — ADAPTIVE kernel\nfilament lands in its own cells only", fontsize=11)

# (c) receptor
ax = axes[2]
ax.plot(fil_x, fil_y, ".", color="k", ms=3, zorder=4)
for x, y in zip(fil_x[::6], fil_y[::6]):
    ax.add_patch(Circle((x, y), 0.20, fill=False, ec="#c02020", lw=0.9, zorder=3))
# one stray parcel close to the station → the spike
ax.plot(3.32, 2.18, ".", color="k", ms=5, zorder=4)
ax.add_patch(Circle((3.32, 2.18), 0.13, fill=False, ec="#c02020", lw=1.4, zorder=3))
ax.plot(*stn_xy, "^", color=C_REC, ms=13, mec="k", zorder=5)
ax.annotate("C(x$_s$) = Σ mᵢ·G(x$_s$−xᵢ; hᵢ)\nevaluated AT the station — no grid",
            xy=(2.8, 5.15), ha="center", fontsize=10)
ax.annotate("one lucky/unlucky parcel with\ntiny h decides the value\n→ single-parcel spikes",
            xy=(3.7, 0.55), ha="center", fontsize=10, color="#c02020")
ax.set_title("(c) RECEPTORS — gridless point KDE\nhonest but high-variance near the source",
             fontsize=11)

fig.suptitle("The four estimator cells — drawn in the geometry of ETEX station 8021: plume core ~20 km (≈2 cells) north of the station",
             fontsize=12, y=1.00)
fig.tight_layout(rect=[0, 0, 1, 0.95])
fig.savefig(os.path.join(OUT, "fig1_estimators_schematic.png"), dpi=150)
plt.close(fig)

# ═══ FIG 2 — the 2×2 RANK decomposition, ETEX vs CAPTEX ═══════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=True)
TERMS = [("R²", lambda s: s["r"] ** 2, "#0072B2"),
         ("1−|FB|/2", lambda s: 1 - abs(s["fb"]) / 2, "#56B4E9"),
         ("FMS/100", lambda s: s["fms"] / 100, "#009E73"),
         ("1−KSP/100", lambda s: 1 - s["ksp"] / 100, "#CC79A7")]
for ax, stats, tag, nrows in ((axes[0], stats_e, "ETEX-1 (400k, 72 h)", len(etex)),
                              (axes[1], stats_c, "CAPTEX R1 (400k, 60 h)", len(captex))):
    xs = np.arange(4)
    bottom = np.zeros(4)
    for tname, tf, tc in TERMS:
        vals = np.array([tf(stats[k]) for k, _, _ in CELLS])
        ax.bar(xs, vals, 0.62, bottom=bottom, color=tc,
               label=tname if ax is axes[0] else None)
        bottom += vals
    for x, (k, _, _) in zip(xs, CELLS):
        s = stats[k]
        ax.text(x, bottom[x] + 0.06, f"RANK {s['rank']:.2f}", ha="center", fontsize=10.5,
                fontweight="bold")
    ax.set_xticks(xs)
    ax.set_xticklabels([lbl.replace(" (", "\n(") +
                        f"\nR={stats[k]['r']:.2f}  FB={stats[k]['fb']:+.2f}"
                        for k, lbl, _ in CELLS], fontsize=8.5)
    ax.set_ylim(0, 3.3)
    ax.set_title(f"{tag} — {nrows} station-windows", fontsize=11)
    ax.grid(axis="y", alpha=0.25)
axes[0].set_ylabel("DATEM RANK = R² + (1−|FB|/2)\n+ FMS/100 + (1−KSP/100)")
fig.legend(loc="upper center", ncol=4, fontsize=9.5, bbox_to_anchor=(0.5, 0.92),
           framealpha=0.95)
fig.suptitle("The 2×2: the cells differ almost only through R² (dark blue) and FB (light blue) — FMS and KSP are flat.\n"
             "On ETEX both moving terms trace to the SAME few near-source pairs (fig. 3–4).",
             fontsize=11)
fig.tight_layout(rect=[0, 0.02, 1, 0.86])
fig.savefig(os.path.join(OUT, "fig2_rank_2x2.png"), dpi=150)
plt.close(fig)

# ═══ FIG 3 — pair-level attribution: it is ONE near-source station ════════════
dist_e = col(etex, "dist")
obs_v = col(etex, "obs")
gu, ga = col(etex, "grid_uni"), col(etex, "grid_ab")
band_idx = np.digitize(dist_e, BAND_EDGES) - 1          # 0..3 inside, 4 beyond

fig = plt.figure(figsize=(13.5, 4.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 1.0])

# (a) uniform vs adaptive per pair, ETEX
ax = fig.add_subplot(gs[0])
eps = 0.3
for b in range(4):
    m = band_idx == b
    ax.loglog(ga[m] + eps, gu[m] + eps, ".", ms=5, color=BAND_COLS[b],
              label=f"{BAND_EDGES[b]}–{BAND_EDGES[b+1]} km (n={m.sum()})", alpha=0.6)
m8021 = col(etex, "stn") == 8021
ax.loglog(ga[m8021] + eps, gu[m8021] + eps, "o", ms=9, mfc="none", mec="k", mew=1.4,
          label="station 8021 windows")
lim = [eps, 1.2e5]
ax.plot(lim, lim, "k-", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("grid · ADAPTIVE kernel (pg/m³)")
ax.set_ylabel("grid · UNIFORM kernel (pg/m³)")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("(a) ETEX: same run pair-by-pair —\nkernels agree except near the source", fontsize=10)

# (b) band means, ETEX
ax = fig.add_subplot(gs[1])
xs = np.arange(4); w = 0.26
means_o = [obs_v[band_idx == b].mean() for b in range(4)]
means_u = [gu[band_idx == b].mean() for b in range(4)]
means_a = [ga[band_idx == b].mean() for b in range(4)]
ax.bar(xs - w, means_o, w, color="0.35", label="observed")
ax.bar(xs, means_u, w, color=C_UNI, label="grid · uniform")
ax.bar(xs + w, means_a, w, color=C_AB, label="grid · adaptive")
ax.set_yscale("log"); ax.set_ylim(0.3, 3e5)
for b in range(4):
    m = band_idx == b
    r_u = np.corrcoef(obs_v[m], gu[m])[0, 1]
    r_a = np.corrcoef(obs_v[m], ga[m])[0, 1]
    ax.text(b, 4e3, f"R {r_u:.2f} / {r_a:.2f}", ha="center", fontsize=8.5)
    if means_a[b] < 2:
        ax.text(b + w, means_a[b] * 1.4, f"{means_a[b]:.1f}", ha="center", fontsize=8, color=C_AB)
ax.set_xticks(xs)
ax.set_xticklabels([f"{BAND_EDGES[b]}–{BAND_EDGES[b+1]} km\nn={np.sum(band_idx==b)}"
                    for b in range(4)], fontsize=9)
ax.set_ylabel("band-mean concentration (pg/m³)")
ax.legend(fontsize=9, loc="upper right", framealpha=0.95)
ax.set_title("(b) ETEX band means (per-band R: uniform / adaptive)\nonly 0–150 km differs — and only in the mean", fontsize=10)

# (c) same scatter for CAPTEX — the control
ax = fig.add_subplot(gs[2])
dist_c = col(captex, "dist")
guc, gac = col(captex, "grid_uni"), col(captex, "grid_ab")
bc = np.digitize(dist_c, BAND_EDGES) - 1
for b in range(4):
    m = bc == b
    if m.sum():
        ax.loglog(gac[m] + eps, guc[m] + eps, ".", ms=5, color=BAND_COLS[b], alpha=0.6,
                  label=f"{BAND_EDGES[b]}–{BAND_EDGES[b+1]} km (n={m.sum()})")
lim = [eps, 1.2e5]
ax.plot(lim, lim, "k-", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("grid · ADAPTIVE kernel (pg/m³)")
ax.set_ylabel("grid · UNIFORM kernel (pg/m³)")
ax.legend(fontsize=8, loc="upper left")
ax.set_title("(c) CAPTEX R1 control: nearest station ≳100 km\n→ kernels indistinguishable everywhere", fontsize=10)

fig.suptitle("Pair-level attribution: the entire ETEX 2×2 contrast lives in the 0–150 km band (74 pairs, dominated by station 8021)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, "fig3_attribution_bands.png"), dpi=150)
plt.close(fig)

# ═══ FIG 4 — the mechanism at station 8021 (real fields) ══════════════════════
# station coordinates + its biggest observed window
r8021 = [r for r in etex if r["stn"] == 8021]
r8021.sort(key=lambda r: -r["obs"])
stn_lat, stn_lon = r8021[0]["lat"], r8021[0]["lon"]
meta = r8021[0]
ti = g_uni_e.frame_of[meta["stop"]]

fig, axes = plt.subplots(1, 3, figsize=(14.0, 5.2), constrained_layout=True,
                         gridspec_kw=dict(width_ratios=[1, 1, 1.2]))
lo0, lo1 = stn_lon - 1.0, stn_lon + 1.0
la0, la1 = stn_lat - 0.55, stn_lat + 0.85
i0, i1 = np.searchsorted(g_uni_e.lon, [lo0, lo1])
j0, j1 = np.searchsorted(g_uni_e.lat, [la0, la1])
lon_e = np.append(g_uni_e.lon[i0:i1] - g_uni_e.dlon / 2, g_uni_e.lon[i1 - 1] + g_uni_e.dlon / 2)
lat_e = np.append(g_uni_e.lat[j0:j1] - g_uni_e.dlat / 2, g_uni_e.lat[j1 - 1] + g_uni_e.dlat / 2)
vmax = max(g.conc[ti, g.klay, j0:j1, i0:i1].max() * 1e15 for g in (g_uni_e, g_ab_e))
norm = LogNorm(vmin=10.0, vmax=vmax)

for cx, (g, gname, gc) in enumerate(((g_uni_e, "UNIFORM kernel", C_UNI),
                                     (g_ab_e, "ADAPTIVE kernel", C_AB))):
    ax = axes[cx]
    f = g.conc[ti, g.klay, j0:j1, i0:i1] * 1e15
    pm = ax.pcolormesh(lon_e, lat_e, np.where(f > 0, f, np.nan), norm=norm,
                       cmap="magma_r", shading="flat")
    for x in lon_e:                                     # 0.1° cell edges
        ax.plot([x, x], [la0, la1], color="0.88", lw=0.3, zorder=2)
    for y in lat_e:
        ax.plot([lo0, lo1], [y, y], color="0.88", lw=0.3, zorder=2)
    ax.plot(ETEX_SRC[1], ETEX_SRC[0], "k*", ms=15, zorder=5, label="source (Monterfil)")
    ax.plot(stn_lon, stn_lat, "^", color="#009E73", ms=12, mec="k", zorder=5,
            label="station 8021")
    ax.axvline(stn_lon, color="0.3", lw=0.8, ls=":", zorder=4)
    mod = g.sample_pg(stn_lon, stn_lat, ti)
    ax.set_title(f"({'ab'[cx]}) {gname}\nmodel at station {mod:,.0f} pg/m³  (obs {meta['obs']:,.0f})",
                 fontsize=10, color=gc)
    ax.set_xlim(lo0, lo1); ax.set_ylim(la0, la1)
    ax.set_xlabel("lon (°E)")
    if cx == 0:
        ax.set_ylabel("lat (°N)")
        ax.legend(loc="lower left", fontsize=8.5)
cb = fig.colorbar(pm, ax=axes[:2], shrink=0.85, pad=0.01, fraction=0.05)
cb.set_label("surface conc. (pg/m³)")

# (c) meridional cut through the station longitude
ax = axes[2]
ix = int(round((stn_lon - g_uni_e.lon[0]) / g_uni_e.dlon))
lats = g_uni_e.lat[j0:j1]
peak_lat, peak_v = 0.0, 0.0
for g, gname, gc in ((g_uni_e, "uniform kernel", C_UNI), (g_ab_e, "adaptive kernel", C_AB)):
    prof = g.conc[ti, g.klay, j0:j1, ix] * 1e15
    if prof.max() > peak_v:
        peak_v, peak_lat = prof.max(), lats[int(prof.argmax())]
    ax.semilogy(lats, np.maximum(prof, 0.05), drawstyle="steps-mid", color=gc, lw=2,
                label=gname)
ax.axvline(stn_lat, color="#009E73", lw=1.6)
ax.plot(stn_lat, meta["obs"], "^", color="#009E73", ms=12, mec="k", zorder=6,
        label=f"station 8021: obs {meta['obs']:,.0f}")
# what the grid estimator actually reports at the station (bilinear sample)
for g, gc, dy in ((g_uni_e, C_UNI, 3.0), (g_ab_e, C_AB, 0.33)):
    mod = g.sample_pg(stn_lon, stn_lat, ti)
    ax.plot(stn_lat, mod, "o", color=gc, ms=8, mec="k", zorder=6)
    ax.text(stn_lat - 0.04, mod * dy, f"grid reads {mod:,.0f}", color=gc,
            fontsize=9, ha="right", va="center")
ax.annotate("", xy=(peak_lat, 2.5e3), xytext=(stn_lat, 2.5e3),
            arrowprops=dict(arrowstyle="<->", lw=1.2, color="0.2"))
ax.text((peak_lat + stn_lat) / 2, 5.5e3, "≈20 km", ha="center", fontsize=9, color="0.2")
ax.set_xlim(la0, la1)
ax.set_ylim(0.05, 4e5)
ax.set_xlabel("lat (°N) along the dotted line in (a)/(b)")
ax.set_ylabel("surface conc. (pg/m³)")
ax.legend(fontsize=8.5, loc="center right")
ax.set_title(f"(c) south–north cut at the station longitude\nframe ending {meta['stop']:%d %b %H:%M}", fontsize=10)

fig.suptitle("The 8021 mechanism (same run, biggest observed window): the filament passes ~0.2° north of the station —\n"
             "the uniform kernel's ±1-cell smearing bridges the gap (false 'hit'), the adaptive kernel reports the miss",
             fontsize=11)
fig.savefig(os.path.join(OUT, "fig4_8021_maps.png"), dpi=150)
plt.close(fig)

# ═══ FIG 5 — the receptor: variance near the source + the N-sweep ═════════════
fig = plt.figure(figsize=(13.5, 4.8))
gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.0, 1.0])

# (a) receptor obs-vs-model, both runs, 8021 highlighted
ax = fig.add_subplot(gs[0])
ru, ra = col(etex, "rec_uni"), col(etex, "rec_ab")
eps = 0.3
ax.loglog(obs_v + eps, ru + eps, ".", ms=4, color=C_REC, alpha=0.45,
          label=f"receptors, uniform run (R={stats_e['rec_uni']['r']:.2f})")
ax.loglog(obs_v + eps, ra + eps, ".", ms=4, color="#c02020", alpha=0.45,
          label=f"receptors, adaptive run (R={stats_e['rec_ab']['r']:.2f})")
ax.loglog(obs_v[m8021] + eps, ru[m8021] + eps, "o", ms=10, mfc="none", mec=C_REC, mew=1.8)
ax.loglog(obs_v[m8021] + eps, ra[m8021] + eps, "o", ms=10, mfc="none", mec="#c02020", mew=1.8)
for o_, u_, a_ in zip(obs_v[m8021], ru[m8021], ra[m8021]):
    ax.plot([o_ + eps, o_ + eps], [min(u_, a_) + eps, max(u_, a_) + eps],
            color="0.4", lw=0.9, zorder=1)
lim = [eps, 1e6]
ax.plot(lim, lim, "k-", lw=0.8)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("observed (pg/m³)"); ax.set_ylabel("receptor estimate (pg/m³)")
ax.legend(fontsize=8.5, loc="upper left")
ax.set_title("(a) SAME estimator, two stochastic runs —\ncircles = the 8021 windows that flip pooled R", fontsize=10)

# (b) receptor stats excluding 8021
ax = fig.add_subplot(gs[1])
m_ex = ~m8021
groups = [("all 2629 pairs", ru, ra, obs_v), ("excluding 8021", ru[m_ex], ra[m_ex], obs_v[m_ex])]
xs = np.arange(2); w = 0.3
r_vals_u = [np.corrcoef(o_, u_)[0, 1] for _, u_, a_, o_ in groups]
r_vals_a = [np.corrcoef(o_, a_)[0, 1] for _, u_, a_, o_ in groups]
ax.bar(xs - w / 2, r_vals_u, w, color=C_REC, label="uniform run")
ax.bar(xs + w / 2, r_vals_a, w, color="#c02020", label="adaptive run")
for x, (vu, va) in zip(xs, zip(r_vals_u, r_vals_a)):
    ax.text(x - w / 2, vu + 0.015, f"{vu:.2f}", ha="center", fontsize=9)
    ax.text(x + w / 2, va + 0.015, f"{va:.2f}", ha="center", fontsize=9)
ax.set_xticks(xs); ax.set_xticklabels([g[0] for g in groups])
ax.set_ylabel("receptor Pearson R (with zeros)")
ax.set_ylim(0, 0.75)
ax.legend(fontsize=9)
ax.set_title("(b) drop ONE station → the two runs agree:\nthe R gap is 8021 point-sampling luck", fontsize=10)

# (c) N-sweep: receptor native resolution + FB converge; grid R flat
ax = fig.add_subplot(gs[2])
sweep = []
for N, gfile, rfile in [(50, "grid_etex_n50k.nc", "receptors_etex_n50k.csv"),
                        (100, "grid_etex_n100k.nc", "receptors_etex_n100k.csv"),
                        (200, "grid_etex_n200k.nc", "receptors_etex_n200k.csv"),
                        (400, "grid_output_etex_ab.nc", "receptors_output_etex_ab.csv"),
                        (800, "grid_etex_n800k.nc", "receptors_etex_n800k.csv")]:
    g = Grid(os.path.join(ETEX_DIR, gfile))
    rcp = Receptors(os.path.join(ETEX_DIR, rfile), g.t0, g.avg_s)
    O, Mg_, Mr_ = [], [], []
    for ob in obs_e:
        vg, vr = g.window_pg(ob), rcp.window_pg(ob)
        if np.isfinite(vg) and np.isfinite(vr):
            O.append(ob["conc"]); Mg_.append(vg); Mr_.append(vr)
    sg, sr = datem_stats(O, Mg_), datem_stats(O, Mr_)
    res = np.array([v for v in rcp.res.values() if np.isfinite(v)])
    sweep.append((N, sg["r"], sr["fb"], np.median(res) / 1e3))
    print(f"  N-sweep {N}k: grid R={sg['r']:.3f} RANK={sg['rank']:.2f} | "
          f"recept FB={sr['fb']:+.2f} R(x_s) p50={np.median(res)/1e3:.1f} km")
Ns = np.array([s[0] for s in sweep], float)
ax.plot(Ns, [s[3] for s in sweep], "o-", color=C_REC, label="receptor native res. p50 (km)")
ax.plot(Ns, 11.7 * (Ns / 50) ** -0.5, "--", color="0.5", lw=1, label="N$^{-1/2}$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xticks(Ns); ax.set_xticklabels([f"{int(n)}k" for n in Ns])
ax.minorticks_off()
ax.set_xlabel("particles N"); ax.set_ylabel("km", color=C_REC)
ax2 = ax.twinx()
ax2.plot(Ns, [s[1] for s in sweep], "s-", color=C_AB, label="grid R (adaptive)")
ax2.plot(Ns, [abs(s[2]) for s in sweep], "^-", color="#7b5cc0", label="receptor |FB|")
ax2.set_ylim(0, 1.0)
ax2.set_ylabel("grid R  /  receptor |FB|")
h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="lower left")
ax.set_title("(c) N-sweep: grid stats converged at 50k;\nreceptor sharpens as N$^{-1/2}$, |FB| → 0", fontsize=10)

fig.suptitle("Why the receptor RANK is lower AND swings between runs: point evaluation has high variance exactly where the plume is thin", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(os.path.join(OUT, "fig5_receptor_variance_nsweep.png"), dpi=150)
plt.close(fig)

print(f"\nfigures + joined pair CSVs → {OUT}")
