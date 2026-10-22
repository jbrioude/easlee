#!/usr/bin/env python3
"""
make_synthetic_wrfout.py — a wrfout file with WRF's exact conventions, built from GFS.

    python3 make_synthetic_wrfout.py GFS_f003.grb2 [GFS_f006.grb2 ...] OUTDIR

There is no WRF installation to produce test data for EASLEE's WRF reader, so this
writes one wrfout_d01_<time> per GFS file with the conventions a real WRF run has,
and stores the truth next to it (truth_<time>.npz) so the reader can be checked
level by level:

  grid      Lambert conformal, MAP_PROJ=1, true latitudes 30/60 N, STAND_LON 10 E,
            centre 48 N 10 E, DX = DY = 30 km, 80 x 70 mass points, Earth radius
            6 370 000 m (WRF's), so the map factor runs ~0.95-1.03
  vertical  41 eta interfaces (bottom_top_stag), 40 mass levels, p_top = 50 hPa,
            hydrostatic: p(eta) = eta (ps - ptop) + ptop, geopotential PH+PHB
            integrated on the interfaces from the terrain with the column's Tv
  winds     U on the x-faces (west_east_stag), V on the y-faces (south_north_stag),
            GRID-RELATIVE (rotated by the cone angle at each face), W on the
            interfaces diagnosed as WRF does: eta-dot from mass continuity with
            these winds (zero at the ground and the top), then the geopotential
            equation, so W = V.grad(h) at the ground.  (Until 2026-10-13 W came
            from GFS omega with W = 0 at the ground: not continuity-consistent,
            which the well-mixed test showed as a near-ground excess of the
            DATA, chi 1.79 in 0-100 m after 3 h, not of the reader.)
  thermo    T = theta - 300 (theta with WRF's R/cp = 2/7), P + PB (perturbation +
            base), PH + PHB, QVAPOR = mixing ratio
  2-D       XLAT/XLONG (+ _U/_V), HGT, PSFC, T2, U10/V10 (grid-relative), HFX
            (positive UPWARD, from GFS avg_ishf), UST (from the GFS stress), PBLH
            (a fixed 800 m: WRF's own PBL scheme output is not reproducible here),
            MAPFAC_M/U/V/UX/VY, MUU/MUV/MUT, Times
Needs numpy, netCDF4, xarray + cfgrib (eccodes).
"""
import sys, os, numpy as np, netCDF4 as nc, xarray as xr

R_EARTH = 6370000.0; RD = 287.0; CP = 7.0 * RD / 2.0; RCP = RD / CP; G = 9.81
TL1, TL2, STDLON, CLAT, CLON = 30.0, 60.0, 10.0, 48.0, 10.0
DX = DY = 30000.0; NI, NJ = 80, 70; NKW = 41; PTOP = 5000.0
D2R = np.pi / 180.0

def cone(t1, t2):
    t1r, t2r = abs(t1) * D2R, abs(t2) * D2R
    return (np.log(np.cos(t1r)) - np.log(np.cos(t2r))) / \
           (np.log(np.tan((90 - abs(t1)) * D2R / 2)) - np.log(np.tan((90 - abs(t2)) * D2R / 2)))
N = cone(TL1, TL2)
def rho_of(lat):            # projected distance from the pole [m]
    return R_EARTH * np.cos(TL1 * D2R) / N * (np.tan((90 - lat) * D2R / 2) / np.tan((90 - TL1) * D2R / 2)) ** N
def fwd(lat, lon):
    r = rho_of(lat); a = N * ((lon - STDLON + 180) % 360 - 180) * D2R
    return r * np.sin(a), -r * np.cos(a)
def inv(x, y):
    r = np.hypot(x, y); lon = STDLON + np.arctan2(x, -y) / N / D2R
    chi = 2 * np.arctan((r * N / (R_EARTH * np.cos(TL1 * D2R))) ** (1 / N) * np.tan((90 - TL1) * D2R / 2))
    return 90 - chi / D2R, lon
def mapfac(lat):
    return N * rho_of(lat) / (R_EARTH * np.cos(lat * D2R))
XC, YC = fwd(CLAT, CLON)
def grid_xy(ioff, joff, ni, nj):
    i = np.arange(1, ni + 1) - (NI + 1) / 2 + ioff; j = np.arange(1, nj + 1) - (NJ + 1) / 2 + joff
    return np.meshgrid(XC + i * DX, YC + j * DY)       # (nj, ni)

def load_gfs(f):
    iso = xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath": "",
          "filter_by_keys": {"typeOfLevel": "isobaricInhPa", "stepType": "instant"}})
    sfc = xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath": "",
          "filter_by_keys": {"typeOfLevel": "surface", "stepType": "instant"}})
    avg = xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath": "",
          "filter_by_keys": {"typeOfLevel": "surface", "stepType": "avg"}})
    t2 = xr.open_dataset(f, engine="cfgrib", backend_kwargs={"indexpath": "",
          "filter_by_keys": {"typeOfLevel": "heightAboveGround", "level": 2}})
    vt = np.datetime64(iso["valid_time"].values, "s").astype(object)
    out = dict(p=iso["isobaricInhPa"].values * 100.0, lat=iso["latitude"].values, lon=iso["longitude"].values,
               time=vt)
    for k in ("u", "v", "t", "q", "gh", "w"):
        out[k] = iso[k].values.astype(np.float64)                  # (np, nlat, nlon)
    for k in ("sp", "orog"):
        out[k] = sfc[k].values.astype(np.float64)
    out["shf"] = avg["avg_ishf"].values; out["tau"] = np.hypot(avg["avg_utaua"].values, avg["avg_vtaua"].values)
    out["t2"] = t2["t2m"].values
    return out

def bilin(G_, field, lat, lon):
    """bilinear on GFS's regular grid (lat descending, lon 0..360); field (..., nlat, nlon)"""
    la, lo = G_["lat"], G_["lon"]; dla = la[1] - la[0]; dlo = lo[1] - lo[0]
    fy = (lat - la[0]) / dla; fx = ((lon % 360) - lo[0]) / dlo
    j0 = np.floor(fy).astype(int); i0 = np.floor(fx).astype(int); ty = fy - j0; tx = fx - i0
    i1 = (i0 + 1) % len(lo)
    f = lambda j, i: field[..., j, i]
    return (1 - ty) * ((1 - tx) * f(j0, i0) + tx * f(j0, i1)) + ty * ((1 - tx) * f(j0 + 1, i0) + tx * f(j0 + 1, i1))

ETAW = 1.0 - (np.linspace(0, 1, NKW) ** 1.6)          # 1 → 0, finer near the ground
ETAW[-1] = 0.0
ETAM = 0.5 * (ETAW[:-1] + ETAW[1:])

def column_profiles(G_, lat, lon):
    """GFS profiles at points → mass/interface columns on the eta grid."""
    prof = {k: bilin(G_, G_[k], lat, lon) for k in ("u", "v", "t", "q", "w", "gh")}  # (np, ...)
    ps = bilin(G_, G_["sp"], lat, lon); hgt = bilin(G_, G_["orog"], lat, lon); t2 = bilin(G_, G_["t2"], lat, lon)
    pm = ETAM[:, None] * (ps - PTOP)[None] + PTOP                   # (nk, npt) mass pressures
    pw = ETAW[:, None] * (ps - PTOP)[None] + PTOP
    lp = np.log(G_["p"])[::-1]                                     # ascending ln p
    res = {}
    npt = ps.size
    for k in ("u", "v", "t", "q", "w"):
        a = prof[k].reshape(len(G_["p"]), npt)[::-1]               # ascending p
        o = np.empty((len(ETAM), npt))
        for n in range(npt):
            ok = G_["p"][::-1] < ps.ravel()[n]                     # above-ground GFS levels only
            x = lp[ok]; y = a[ok, n]
            if k == "t":                                           # add the 2 m temperature at the surface
                x = np.append(x, np.log(ps.ravel()[n])); y = np.append(y, t2.ravel()[n])
            o[:, n] = np.interp(np.log(pm[:, n]), x, y)
        res[k] = o
    tv = res["t"] * (1 + 0.608 * res["q"])
    zw = np.empty((NKW, npt)); zw[0] = hgt.ravel()
    for k in range(NKW - 1):
        zw[k + 1] = zw[k] + RD * tv[k] / G * np.log(pw[k] / pw[k + 1])
    res.update(pm=pm, pw=pw, ps=ps.ravel(), hgt=hgt.ravel(), zw=zw, t2=t2.ravel())
    return res

def to_grid_rel(ue, ve, lon):
    a = N * ((lon - STDLON + 180) % 360 - 180) * D2R
    return ue * np.cos(a) - ve * np.sin(a), ue * np.sin(a) + ve * np.cos(a)

def write(G_, outdir):
    tstr = G_["time"].strftime("%Y-%m-%d_%H:%M:%S")
    xm, ym = grid_xy(0, 0, NI, NJ); latm, lonm = inv(xm, ym)
    xu, yu = grid_xy(-0.5, 0, NI + 1, NJ); latu, lonu = inv(xu, yu)
    xv, yv = grid_xy(0, -0.5, NI, NJ + 1); latv, lonv = inv(xv, yv)
    cm = column_profiles(G_, latm.ravel(), lonm.ravel())
    cu = column_profiles(G_, latu.ravel(), lonu.ravel())
    cv = column_profiles(G_, latv.ravel(), lonv.ravel())
    nkm = NKW - 1
    rs = lambda a, nj, ni: a.reshape(a.shape[0], nj, ni)
    # winds: grid-relative at the face points
    ug, _ = to_grid_rel(cu["u"], cu["v"], lonu.ravel()[None]); _, vg = to_grid_rel(cv["u"], cv["v"], lonv.ravel()[None])
    theta = cm["t"] * (1e5 / cm["pm"]) ** RCP
    rhom = cm["pm"] / (RD * cm["t"] * (1 + 0.608 * cm["q"]))
    # W on interfaces as WRF diagnoses it, consistent with these winds:
    #   mass continuity  mu_t + m^2 [d(mu u/m)/dx + d(mu v/m)/dy] + d(mu etadot)/deta = 0,
    #   etadot = 0 at the ground and the top (mu_t closes the column), then the
    #   geopotential equation  g w = m (u dphi/dx + v dphi/dy) + etadot dphi/deta
    #   (dphi/dt left out: two archive times only; EASLEE neglects dz_k/dt too).
    #   At the ground this is w = V.grad(h), WRF's kinematic condition.
    om = cm["w"]
    nkm_ = NKW - 1; deta = ETAW[:-1] - ETAW[1:]                      # > 0, layer k
    mu_m = (cm["ps"] - PTOP).reshape(NJ, NI); m_m = mapfac(latm)
    fu = rs(ug, NJ, NI + 1) * ((cu["ps"] - PTOP) / mapfac(latu).ravel()).reshape(NJ, NI + 1)[None]
    fv = rs(vg, NJ + 1, NI) * ((cv["ps"] - PTOP) / mapfac(latv).ravel()).reshape(NJ + 1, NI)[None]
    div = m_m[None] ** 2 * ((fu[:, :, 1:] - fu[:, :, :-1]) / DX + (fv[:, 1:, :] - fv[:, :-1, :]) / DY)
    mu_t = -(div * deta[:, None, None]).sum(0)
    Om = np.zeros((NKW, NJ, NI))                                     # mu * etadot on the interfaces
    for k in range(nkm_):
        Om[k + 1] = Om[k] + (mu_t + div[k]) * deta[k]
    Om[-1] = 0.0                                                     # (closes to round-off)
    etd = Om / mu_m[None]
    phi = G * rs(cm["zw"], NJ, NI)                                   # (NKW, NJ, NI)
    ugm, vgm = to_grid_rel(cm["u"], cm["v"], lonm.ravel()[None])     # grid-relative, mass points
    ugm = rs(ugm, NJ, NI); vgm = rs(vgm, NJ, NI)
    ui = np.concatenate([ugm[:1], 0.5 * (ugm[1:] + ugm[:-1]), ugm[-1:]])
    vi = np.concatenate([vgm[:1], 0.5 * (vgm[1:] + vgm[:-1]), vgm[-1:]])
    dphidx = np.gradient(phi, DX, axis=2); dphidy = np.gradient(phi, DY, axis=1)
    dphide = np.gradient(phi, ETAW, axis=0)
    W = ((m_m[None] * (ui * dphidx + vi * dphidy) + etd * dphide) / G).reshape(NKW, NJ * NI)
    ue_sfc, ve_sfc = cm["u"][0], cm["v"][0]
    u10, v10 = to_grid_rel(ue_sfc, ve_sfc, lonm.ravel())
    shf = bilin(G_, G_["shf"], latm.ravel(), lonm.ravel()); tau = bilin(G_, G_["tau"], latm.ravel(), lonm.ravel())
    rho_s = cm["ps"] / (RD * cm["t2"])
    ust = np.sqrt(np.maximum(tau, 0) / rho_s)
    path = os.path.join(outdir, f"wrfout_d01_{tstr}")
    d = nc.Dataset(path, "w", format="NETCDF4")
    for name, n in (("Time", None), ("DateStrLen", 19), ("west_east", NI), ("south_north", NJ),
                    ("bottom_top", nkm), ("west_east_stag", NI + 1), ("south_north_stag", NJ + 1),
                    ("bottom_top_stag", NKW)):
        d.createDimension(name, n)
    def var(name, dims, data, stagger="", units=""):
        v = d.createVariable(name, "f4", dims); v.stagger = stagger; v.units = units
        v[:] = data[None].astype(np.float32); return v
    T_ = d.createVariable("Times", "S1", ("Time", "DateStrLen")); T_[0, :] = np.array(list(tstr), dtype="S1")
    M3 = ("Time", "bottom_top", "south_north", "west_east")
    var("U", ("Time", "bottom_top", "south_north", "west_east_stag"), rs(ug, NJ, NI + 1), "X", "m s-1")
    var("V", ("Time", "bottom_top", "south_north_stag", "west_east"), rs(vg, NJ + 1, NI), "Y", "m s-1")
    var("W", ("Time", "bottom_top_stag", "south_north", "west_east"), rs(W, NJ, NI), "Z", "m s-1")
    var("T", M3, rs(theta - 300.0, NJ, NI), "", "K")
    p_tot = rs(cm["pm"], NJ, NI); var("PB", M3, 0.97 * p_tot, "", "Pa"); var("P", M3, 0.03 * p_tot, "", "Pa")
    ph_tot = rs(G * cm["zw"], NJ, NI)
    var("PHB", ("Time", "bottom_top_stag", "south_north", "west_east"), 0.97 * ph_tot, "Z", "m2 s-2")
    var("PH", ("Time", "bottom_top_stag", "south_north", "west_east"), 0.03 * ph_tot, "Z", "m2 s-2")
    var("QVAPOR", M3, rs(cm["q"] / (1 - cm["q"]), NJ, NI), "", "kg kg-1")
    S2 = ("Time", "south_north", "west_east")
    for nm, a in (("XLAT", latm), ("XLONG", lonm)):
        var(nm, S2, a)
    for nm, a in (("XLAT_U", latu), ("XLONG_U", lonu)):
        var(nm, ("Time", "south_north", "west_east_stag"), a, "X")
    for nm, a in (("XLAT_V", latv), ("XLONG_V", lonv)):
        var(nm, ("Time", "south_north_stag", "west_east"), a, "Y")
    r2 = lambda a: a.reshape(NJ, NI)
    var("HGT", S2, r2(cm["hgt"]), "", "m"); var("PSFC", S2, r2(cm["ps"]), "", "Pa"); var("T2", S2, r2(cm["t2"]), "", "K")
    var("U10", S2, r2(u10), "", "m s-1"); var("V10", S2, r2(v10), "", "m s-1")
    var("HFX", S2, r2(shf), "", "W m-2"); var("UST", S2, r2(ust), "", "m s-1")
    var("PBLH", S2, np.full((NJ, NI), 800.0), "", "m")
    mu = r2(cm["ps"] - PTOP)
    var("MUT", S2, mu, "", "Pa"); var("MUU", ("Time", "south_north", "west_east_stag"), (cu["ps"] - PTOP).reshape(NJ, NI + 1), "X", "Pa")
    var("MUV", ("Time", "south_north_stag", "west_east"), (cv["ps"] - PTOP).reshape(NJ + 1, NI), "Y", "Pa")
    var("MAPFAC_M", S2, mapfac(latm)); var("MAPFAC_U", ("Time", "south_north", "west_east_stag"), mapfac(latu), "X")
    var("MAPFAC_V", ("Time", "south_north_stag", "west_east"), mapfac(latv), "Y")
    var("MAPFAC_UX", ("Time", "south_north", "west_east_stag"), mapfac(latu), "X")
    var("MAPFAC_VY", ("Time", "south_north_stag", "west_east"), mapfac(latv), "Y")
    for k, v in dict(DX=DX, DY=DY, CEN_LAT=CLAT, CEN_LON=CLON, TRUELAT1=TL1, TRUELAT2=TL2, STAND_LON=STDLON,
                     MOAD_CEN_LAT=CLAT, MAP_PROJ=1, GRID_ID=1, PARENT_ID=0, I_PARENT_START=1, J_PARENT_START=1,
                     **{"WEST-EAST_PATCH_END_UNSTAG": NI, "SOUTH-NORTH_PATCH_END_UNSTAG": NJ,
                        "BOTTOM-TOP_PATCH_END_STAG": NKW, "WEST-EAST_GRID_DIMENSION": NI + 1,
                        "SOUTH-NORTH_GRID_DIMENSION": NJ + 1, "BOTTOM-TOP_GRID_DIMENSION": NKW}).items():
        d.setncattr(k, np.float32(v) if isinstance(v, float) else np.int32(v))
    d.TITLE = "SYNTHETIC wrfout from GFS (validation/make_synthetic_wrfout.py) - not a WRF run"
    d.close()
    # truth at mass points: physical T [K], p [Pa], mass-level height AGL [m], EARTH-relative u, v, q (specific)
    zm = 0.5 * (cm["zw"][1:] + cm["zw"][:-1]) - cm["hgt"][None]
    tdir = os.path.join(os.path.dirname(os.path.abspath(outdir)), os.path.basename(os.path.abspath(outdir)) + "_truth")
    os.makedirs(tdir, exist_ok=True)                   # NOT in the met directory (the scan would open it)
    t = nc.Dataset(os.path.join(tdir, f"truth_{tstr}.nc"), "w")
    t.createDimension("k", nkm); t.createDimension("j", NJ); t.createDimension("i", NI)
    t.createDimension("kw", NKW)                       # interfaces (w-levels)
    for name, a in dict(zw=rs(cm["zw"] - cm["hgt"][None], NJ, NI), wi=rs(W, NJ, NI)).items():
        t.createVariable(name, "f8", ("kw", "j", "i"))[:] = a          # interface height AGL, interface W
    for name, a in dict(T=rs(cm["t"], NJ, NI), p=p_tot, zm=rs(zm, NJ, NI), u=rs(cm["u"], NJ, NI),
                        v=rs(cm["v"], NJ, NI), q=rs(cm["q"], NJ, NI), w=rs(0.5 * (W[1:] + W[:-1]), NJ, NI)).items():
        t.createVariable(name, "f8", ("k", "j", "i"))[:] = a          # EARTH-relative u, v; zm above ground
    for name, a in dict(lat=latm, lon=lonm, mapfac=mapfac(latm), hgt=cm["hgt"].reshape(NJ, NI)).items():
        t.createVariable(name, "f8", ("j", "i"))[:] = a
    t.cone = N; t.stand_lon = STDLON
    t.close()
    print("wrote", path, "map factor", mapfac(latm).min().round(4), "-", mapfac(latm).max().round(4))

if __name__ == "__main__":
    *files, outdir = sys.argv[1:]
    os.makedirs(outdir, exist_ok=True)
    for f in files:
        write(load_gfs(f), outdir)
