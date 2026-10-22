#!/usr/bin/env python3
# check_etadot_levels.py FILE.grb [...] — where does the archived eta-dot sit?
#
# Mass continuity in hybrid coordinates fixes the flux through every interface
# from u, v and the surface pressure alone (no time derivative):
#     F(k+1/2) = eta-dot dp/deta at the interface below layer k (top-down)
#              = B(k+1/2) * sum_all D_j - sum_{j<=k} D_j,   D_j = div(v_j dp_j)
# The archived value at GRIB level k is correlated with F at the interface
# ABOVE (k-1/2), BELOW (k+1/2) and with their mean (the full level k).  The
# right placement gives the same slope at every level; the wrong ones drift
# toward the ground.  Units are detected from the magnitude (eta-dot [s^-1]
# ~1e-6 is multiplied by dp/deta of its layer; flex_extract eta-dot*dp/deta
# [Pa/s] ~1e-1 is used as is).  Global (wrapping) or regional regular lat-lon.
#
# Measured 2026-10-22: flex_extract (EAg, ERA5 global) → BELOW (half levels,
# level 137 = the ground); raw MARS (ETEX IFS_EA) → MEAN (full levels).
# The reader decides the same thing from rms(137)/rms(136) (_etadot_levels!).
#
# Needs eccodes + numpy (on zion: PYTHONPATH=/data/jbrioude/met_samples/pysite).
import sys
import numpy as np
import eccodes as ec

def load(fn):
    U, V, E = {}, {}, {}
    sp = lnsp = pv = None; meta = None
    with open(fn, "rb") as f:
        while (g := ec.codes_grib_new_from_file(f)) is not None:
            n = ec.codes_get(g, "shortName"); tl = ec.codes_get(g, "typeOfLevel")
            Ni = ec.codes_get(g, "Ni"); Nj = ec.codes_get(g, "Nj")
            if tl == "hybrid" and n in ("u", "v", "etadot"):
                a = ec.codes_get_values(g).reshape(Nj, Ni)
                {"u": U, "v": V, "etadot": E}[n][ec.codes_get(g, "level")] = a
                if pv is None:
                    pv = ec.codes_get_array(g, "pv")
                    meta = (Ni, Nj, ec.codes_get(g, "latitudeOfFirstGridPointInDegrees"),
                            ec.codes_get(g, "iDirectionIncrementInDegrees"),
                            ec.codes_get(g, "jDirectionIncrementInDegrees"))
            elif n == "sp":
                sp = ec.codes_get_values(g).reshape(Nj, Ni)
            elif n == "lnsp":
                lnsp = ec.codes_get_values(g).reshape(Nj, Ni)
            ec.codes_release(g)
    if sp is None:
        sp = np.exp(lnsp)
    return U, V, E, sp, pv, meta

def main(fn):
    U, V, E, sp, pv, (Ni, Nj, la0, di, dj) = load(fn)
    L = sorted(U); nl = len(L)
    A = pv[:len(pv) // 2]; B = pv[len(pv) // 2:]; p0 = 101325.0
    glob = abs(Ni * di - 360.0) < 1e-6
    R = 6.371e6; dl = np.deg2rad(di); dphi = np.deg2rad(dj)
    lat = np.deg2rad(la0 - dj * np.arange(Nj))[:, None]; cl = np.maximum(np.cos(lat), 1e-3)
    def div(fu, fv):
        out = np.full(fu.shape, np.nan)
        if glob:
            dx = (np.roll(fu, -1, 1) - np.roll(fu, 1, 1)) / (2 * dl * R * cl)
        else:
            dx = np.full(fu.shape, np.nan)
            dx[:, 1:-1] = (fu[:, 2:] - fu[:, :-2]) / (2 * dl * R * cl)
        fvc = fv * cl
        out[1:-1] = dx[1:-1] + (fvc[:-2] - fvc[2:]) / (2 * dphi * R * cl[1:-1])   # rows run north → south
        return out
    D = np.zeros((nl, Nj, Ni))
    for i, k in enumerate(L):
        dp = (A[k] - A[k - 1]) + (B[k] - B[k - 1]) * sp
        D[i] = div(U[k] * dp, V[k] * dp)
    tot = D.sum(0); cum = np.cumsum(D, 0)
    F = np.zeros((nl + 1, Nj, Ni))
    for h in range(1, nl + 1):
        F[h] = B[h] * tot - cum[h - 1]
    mask = np.isfinite(tot) & (np.abs(np.rad2deg(lat)) < 60.0)
    if not glob:
        mask[:3] = mask[-3:] = False; mask[:, :3] = mask[:, -3:] = False
    per_s = np.sqrt(np.mean(E[L[nl // 2]][mask] ** 2)) < 1e-3
    print(f"{fn}: {'eta-dot [s^-1] x dp/deta' if per_s else 'eta-dot*dp/deta [Pa/s]'}, {'global' if glob else 'regional'}")
    def st(x, y):
        x = x[mask]; y = y[mask]
        return np.corrcoef(x, y)[0, 1], np.dot(x, y) / np.dot(x, x)
    print(f"{'level':>5s} {'corr above':>10s} {'slope':>6s} {'corr below':>10s} {'slope':>6s} {'corr mean':>9s} {'slope':>6s}  rms [Pa/s]")
    for k in [k for k in (137, 136, 135, 133, 130, 125, 120, 110, 100, 90, 70) if k in E]:
        e = E[k]
        if per_s:
            e = e * ((A[k] - A[k - 1]) + (B[k] - B[k - 1]) * sp) / ((A[k] - A[k - 1]) / p0 + (B[k] - B[k - 1]))
        ca, sa = st(F[k - 1], e)
        cb, sb = st(F[k], e) if k < L[-1] else (np.nan, np.nan)
        cm, sm = st(0.5 * (F[k - 1] + F[k]), e)
        print(f"{k:5d} {ca:10.3f} {sa:6.2f} {cb:10.3f} {sb:6.2f} {cm:9.3f} {sm:6.2f}  {np.sqrt(np.mean(e[mask] ** 2)):.3e}")
    r = np.sqrt(np.mean(E[L[-1]][mask] ** 2)) / np.sqrt(np.mean(E[L[-2]][mask] ** 2))
    print(f"rms(level {L[-1]}) / rms(level {L[-2]}) = {r:.3f}  (< 0.15: half levels; the reader's auto test)")

if __name__ == "__main__":
    with np.errstate(invalid="ignore", divide="ignore"):
        for fn in sys.argv[1:]:
            main(fn)
