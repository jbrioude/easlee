# Per-model derived-field calculations for the multi-source met readers.
# Each function brings the raw loaded fields to the common EASLEE convention:
#   k = 1 surface slot (zlev = 0), k = 2..Nk ascending levels;
#   rho, theta, drhodz, levdot [level-index/s], puv [Pa] filled;
#   tke/tls from Hanna when surface parameters are available, zero otherwise.
#
# Vertical mixing (Fix C, ADAPTIVE_KERNEL §8): the "tke" field's MEANING is a
# per-driver property published by the reader —
#   Hanna-fill drivers here (MERRA-2 / GFS / AROME-HDF, like the ERA5 GRIB
#   path): tke := 0.5·σ_w,H², so σ_w = √(2·tke) and K_z = σ_w²·τ_Lw,H;
#   MesoNH (and AROME native e): tke := the archived full TKE e; σ_w² = 2·e·r_w
#   and K_z = C_K·ℓ·√e (CBR), constants from the listvariables namelist.
# The closure tuple is derived in gridcheck_onefile! (vmix_namelist.jl) and
# applied in physics/turbulence.jl (tke_partition) — nothing here branches on
# it, each reader just keeps producing its declared field.
#
# Loader conventions on entry (see _read_netcdf_file! / _read_grib_messages!):
#   MERRA2  : mid-level vars at k=1..Nk-1 bottom-up (flipped by the loader);
#             edge vars (CMFMC, PF*) at k=1..Nk bottom-up (edge 1 = surface).
#   GFS     : isobaric-ladder vars at k=1..Nk-1 bottom-up (1000 hPa first).
#   AROME_HDF: height-ladder vars at k=1..Nk-1 bottom-up (20 m first).
#   MESONH  : native MesoNH levels at k=1..Nk bottom-up (k=1 is the below-
#             ground ghost level, which becomes the surface slot).

const _RD  = 287.05f0
const _G0  = 9.81f0

# Shift mid-level 3D fields up one slot to free k=1 for the surface layer.
function _shift_up_one!(metfields3D, keys, Nk, indt)
   for key in keys
      haskey(metfields3D, key) || continue
      metfields3D[key][:,:,2:Nk,indt] = metfields3D[key][:,:,1:Nk-1,indt]
   end
   return nothing
end

# rho (moist), theta, drhodz from puv/T/qv already on the common layout.
function _thermo_from_pTq!(metfields3D, Nk, indt)
   p  = metfields3D["puv"]
   T  = metfields3D["T"]
   qv = haskey(metfields3D, "qv") ? metfields3D["qv"] : nothing
   tv = qv === nothing ? T[:,:,:,indt] :
        T[:,:,:,indt] .* (1.0f0 .+ 0.608f0 .* max.(qv[:,:,:,indt], 0.0f0))
   metfields3D["rho"][:,:,:,indt] = p[:,:,:,indt] ./ (_RD .* max.(tv, 100.0f0))
   metfields3D["theta"][:,:,:,indt] = T[:,:,:,indt] .*
        (100000.0f0 ./ max.(p[:,:,:,indt], 1.0f0)).^0.28573f0
   rho = metfields3D["rho"]; zl = metfields3D["zlev"]
   for k in 2:Nk-1
      metfields3D["drhodz"][:,:,k,indt] = (rho[:,:,k+1,indt] .- rho[:,:,k-1,indt]) ./
          max.(zl[:,:,k+1,indt] .- zl[:,:,k-1,indt], 1.0f0)
   end
   metfields3D["drhodz"][:,:,1,indt]  = metfields3D["drhodz"][:,:,2,indt]
   metfields3D["drhodz"][:,:,Nk,indt] = metfields3D["drhodz"][:,:,Nk-1,indt]
   fill!(view(metfields3D["dsigwdz"], :, :, :, indt), 0.0f0)
   return nothing
end

# Surface turbulence parameters → Hanna tke/tls.
# Derivation chain (first available source wins):
#   ustar : archived → surface stress (tauu/tauv [N/m²]) → neutral log-law
#           from the lowest-level wind (z0 from lsm when present)
#   ol    : archived → sensible heat flux shf [W/m², +up] → surface-layer
#           bulk Richardson number
#   pblh  : archived → bulk-Richardson profile (Troen & Mahrt, Ri_c = 0.25 —
#           the same diagnostic ERA5 uses for blh)
# Falls back to zeroed tke/tls (σw floors in the turbulence scheme) only when
# even the lowest-level winds are unusable.
# Obukhov length from the surface-layer bulk Richardson number, STACKED layout:
# θ and u, v on the full levels of slots 2 and 3, at heights zf(s) = ½(zlev[s-1] +
# zlev[s]); Rib = g Δθ Δz / (θ ΔU²) between them, ζ = Rib (unstable) or
# Rib/(1 − 5 Rib) (stable) at the upper of the two, L = zf(3)/ζ.  Used when no
# heat flux is archived (2026-10-22: the ERA5 path then left L = 0 — neutral
# turbulence everywhere — on every ERA5_CAPTEX/CONV/FULL/CAPE run).
function _ol_bulk_richardson!(metfields3D, metfields2D, indt, Ni, Nj)
   th = metfields3D["theta"]; u = metfields3D["Ux"]; v = metfields3D["Vy"]
   zl = metfields3D["zlev"]; ol = metfields2D["ol"]
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         z2 = 0.5f0 * (zl[i,j,1,indt] + zl[i,j,2,indt]); z3 = 0.5f0 * (zl[i,j,2,indt] + zl[i,j,3,indt])
         dz  = z3 - z2
         dth = th[i,j,3,indt] - th[i,j,2,indt]
         du2 = (u[i,j,3,indt] - u[i,j,2,indt])^2 + (v[i,j,3,indt] - v[i,j,2,indt])^2 + 0.01f0
         (isfinite(dth) && isfinite(du2) && dz > 0.5f0) || continue
         rib  = _G0 * dz * dth / (max(th[i,j,2,indt], 200.0f0) * du2)
         zeta = rib >= 0.0f0 ? rib / max(1.0f0 - 5.0f0 * rib, 0.1f0) : rib
         ol[i,j,indt] = abs(zeta) < 1.0f-4 ? 5000.0f0 : clamp(z3 / zeta, -5000.0f0, 5000.0f0)
      end
   end
   return nothing
end

function _surface_turbulence!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   has_ustar, has_pblh = _derive_surface_params!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   if has_ustar && has_pblh
      hanna_fill_tke!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   else
      # Silent degradation was a trap: with zeroed tke/tls the turbulence
      # scheme runs at its σw floor everywhere (≈ laminar transport).
      @warn "surface turbulence: ustar/pblh underivable from this met — tke/tls zeroed, σw floors apply" maxlog=3
      fill!(view(metfields3D["tke"], :, :, :, indt), 0.0f0)
      fill!(view(metfields3D["tls"], :, :, :, indt), 0.0f0)
   end
   return nothing
end

# Derivation of the 2D surface parameters only (no Hanna tke fill) — used
# directly by datasets with native TKE (MesoNH) that still need ustar/ol/pblh
# for the Hanna option-b path and the PBL entrainment zone.
# `native_tke=true` (MesoNH) prefers the TKE-threshold PBL height over the bulk
# Richardson profile — NEVER set it for datasets whose tke array is the Hanna
# 0.5·σw² fill (the 2-slot rotation would feed the previous met time's fill
# back into the diagnostic).
function _derive_surface_params!(metfields3D, metfields2D, indt, Ni, Nj, Nk;
                                 native_tke::Bool=false)
   rho_sfc = metfields3D["rho"][:,:,2,indt]
   T_sfc   = metfields3D["T"][:,:,2,indt]

   _nzmax(a) = maximum(x -> isfinite(x) ? abs(x) : 0.0f0, a)
   has_ustar = haskey(metfields2D, "ustar") && _nzmax(metfields2D["ustar"][:,:,indt]) > 0
   has_pblh  = haskey(metfields2D, "pblh")  && _nzmax(metfields2D["pblh"][:,:,indt])  > 0
   has_ol    = haskey(metfields2D, "ol")    && _nzmax(metfields2D["ol"][:,:,indt]) > 0

   if !has_ustar && haskey(metfields2D, "tauu") && haskey(metfields2D, "tauv")
      tau = sqrt.(metfields2D["tauu"][:,:,indt].^2 .+ metfields2D["tauv"][:,:,indt].^2)
      if _nzmax(tau) > 0
         metfields2D["ustar"][:,:,indt] .= sqrt.(tau ./ max.(rho_sfc, 0.1f0))
         has_ustar = true
      end
   end

   if !has_ustar && haskey(metfields2D, "ustar")
      # Neutral log-law from the lowest model-level wind.
      u  = metfields3D["Ux"]; v = metfields3D["Vy"]; zl = metfields3D["zlev"]
      lsm = haskey(metfields2D, "lsm") ? metfields2D["lsm"] : nothing
      Threads.@threads :static for j in 1:Nj
         @inbounds for i in 1:Ni
            spd = sqrt(u[i,j,2,indt]^2 + v[i,j,2,indt]^2)
            isfinite(spd) || continue
            z1  = max(0.5f0 * zl[i,j,2,indt], 2.0f0)      # stacked: full level 1 at zlev[2]/2
            z0  = lsm === nothing ? 0.01f0 : (lsm[i,j,indt] > 0.5f0 ? 0.05f0 : 2.0f-4)
            metfields2D["ustar"][i,j,indt] = clamp(0.4f0 * spd / log(z1 / z0), 0.02f0, 1.5f0)
         end
      end
      has_ustar = true
   end

   if !has_ol && has_ustar && haskey(metfields2D, "shf") &&
      _nzmax(metfields2D["shf"][:,:,indt]) > 0
      # L = −ρ·cp·T·u*³ / (κ·g·H),  H positive upward
      H  = metfields2D["shf"][:,:,indt]
      u3 = max.(metfields2D["ustar"][:,:,indt], 0.001f0).^3
      ol = -rho_sfc .* 1004.0f0 .* T_sfc .* u3 ./
           (0.4f0 .* _G0 .* (H .+ sign.(H) .* 1.0f-6 .+ 1.0f-6))
      metfields2D["ol"][:,:,indt] .= clamp.(ol, -5000.0f0, 5000.0f0)
      has_ol = true
   end

   if !has_ol && haskey(metfields2D, "ol")
      _ol_bulk_richardson!(metfields3D, metfields2D, indt, Ni, Nj)    # stacked heights
      has_ol = true
   end

   if !has_pblh && haskey(metfields2D, "pblh")
      if native_tke && haskey(metfields3D, "tke") &&
         maximum(x -> isfinite(x) ? x : 0.0f0, view(metfields3D["tke"], :, :, 2, indt)) > 0.01f0
         _pblh_from_tke!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
         @info "pblh: derived from the archived TKE profile (threshold diagnostic)" maxlog=1
      else
         _pblh_from_rib!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
         @info "pblh: derived from the bulk-Richardson profile (Ri_c = 0.25)" maxlog=1
      end
      has_pblh = true
   end

   return has_ustar, has_pblh
end

# Boundary-layer height from a TKE threshold (datasets with native TKE, e.g.
# MesoNH): lowest height where TKE falls below max(5 % of the near-surface
# value, 0.01 m² s⁻²) — the standard prognostic-TKE PBL diagnostic (COSMO /
# HIRLAM convention).  Falls back to the Rib range clamp [50, 6000] m.
function _pblh_from_tke!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   tke = metfields3D["tke"]; zl = metfields3D["zlev"]
   pblh = metfields2D["pblh"]
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         t2 = tke[i,j,2,indt]
         (isfinite(t2) && t2 > 0.01f0) || (pblh[i,j,indt] = 50.0f0; continue)
         thr = max(0.05f0 * t2, 0.01f0)
         h = 50.0f0
         for k in 3:Nk
            tk = tke[i,j,k,indt]
            if !(isfinite(tk) && tk >= thr)
               h = zl[i,j,k,indt]
               break
            end
            zl[i,j,k,indt] > 6000.0f0 && (h = 6000.0f0; break)
         end
         pblh[i,j,indt] = clamp(h, 50.0f0, 6000.0f0)
      end
   end
   return nothing
end

# Boundary-layer height from the bulk Richardson number profile (Troen &
# Mahrt 1986; Ri_c = 0.25 — the same diagnostic the IFS uses for blh):
# lowest height where Rib(z) = g·z·(θv(z)−θv_sfc) / (θv_sfc·|U(z)|²) crosses Ri_c.
# The shear is taken from the SURFACE, where the wind is zero (Vogelezang &
# Holtslag 1996; Seidel et al. 2012; the GFS PBL scheme).  Until 2026-11-02 it
# was U(z) − U(lowest level): on GFS isobaric data (lowest full level ≈ 107 m)
# that drops the shear below it, and pblh came out 4× below GFS's own hpbl over
# the sea (medians 175 vs 734 m; 722 m with the surface reference).
function _pblh_from_rib!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   th = metfields3D["theta"]; u = metfields3D["Ux"]; v = metfields3D["Vy"]
   zl = metfields3D["zlev"]
   qv = haskey(metfields3D, "qv") ? metfields3D["qv"] : nothing
   pblh = metfields2D["pblh"]
   Ric  = 0.25f0
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         thvs = th[i,j,2,indt] * (qv === nothing ? 1.0f0 : 1.0f0 + 0.61f0 * max(qv[i,j,2,indt], 0.0f0))
         isfinite(thvs) || (pblh[i,j,indt] = 200.0f0; continue)
         h = 100.0f0
         # stacked layout: θ, u, v of slot k are full-level values at ½(zlev[k-1] + zlev[k])
         ribp = 0.0f0; zp = 0.5f0 * (zl[i,j,1,indt] + zl[i,j,2,indt])
         for k in 3:Nk
            z = 0.5f0 * (zl[i,j,k-1,indt] + zl[i,j,k,indt])
            thv = th[i,j,k,indt] * (qv === nothing ? 1.0f0 : 1.0f0 + 0.61f0 * max(qv[i,j,k,indt], 0.0f0))
            du2 = u[i,j,k,indt]^2 + v[i,j,k,indt]^2 + 0.1f0
            rib = _G0 * z * (thv - thvs) / (max(thvs, 200.0f0) * du2)
            if rib >= Ric
               # interpolate the crossing between the previous level and here
               f = (Ric - ribp) / max(rib - ribp, 1.0f-6)
               h = zp + f * (z - zp)
               break
            end
            ribp = rib; zp = z
            z > 6000.0f0 && (h = 6000.0f0; break)
         end
         pblh[i,j,indt] = clamp(h, 50.0f0, 6000.0f0)
      end
   end
   return nothing
end

# Pseudo-3-D precipitation flux from a surface precipitation RATE [kg/m²/s]:
# the rate is applied at every level below the cloud top (highest level with
# cloud fraction > 0.2), so in-cloud/below-cloud wet scavenging sees rain
# falling through the whole column under the cloud.  Used for datasets that
# archive only surface precipitation (GFS prate).
function _tpf_from_surface_precip!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   (haskey(metfields3D, "tpf") && haskey(metfields2D, "prate")) || return nothing
   tpf = metfields3D["tpf"]; zl = metfields3D["zlev"]
   cc  = haskey(metfields3D, "cc") ? metfields3D["cc"] : nothing
   pr  = metfields2D["prate"]
   fill!(view(tpf, :, :, :, indt), 0.0f0)
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         p = pr[i,j,indt]
         (isfinite(p) && p > 1.0f-7) || continue
         ztop = 5000.0f0
         if cc !== nothing
            for k in Nk:-1:2
               zl[i,j,k,indt] > 14000.0f0 && continue
               if cc[i,j,k,indt] > 0.2f0
                  ztop = max(zl[i,j,k,indt], 1000.0f0)
                  break
               end
            end
         end
         for k in 1:Nk
            zl[i,j,k,indt] <= ztop || break
            tpf[i,j,k,indt] = p
         end
      end
   end
   return nothing
end

# ── Vertical velocity on terrain-following levels ────────────────────────────
# The advance keeps a parcel's level index while it moves horizontally, so it
# already follows the model level over the terrain; levdot must therefore be the
# COORDINATE velocity
#     dk/dt = (w − u·∂z_k/∂x − v·∂z_k/∂y) / (∂z/∂k),     z_k = zlev_k + z_s (ASL),
# not the physical w (FLEXPART-WRF makes the same correction).  Using w counted
# the terrain-forced part twice: near the ground over a slope w ≈ V·∇z_s, which
# the level already provides (2026-10-13; La Réunion AROME: 0.25 m/s mean over
# land).  The tendency ∂z_k/∂t of the (pressure-based) levels is neglected.
# u, v: mass-point winds along the grid axes.  Slopes over the grid spacing
# measured on the sphere from the nest's own lon/lat (any projection, map factor
# included), cached per nest.
function _grid_spans_m!(metconfig, Ni, Nj)
   c = get(metconfig, "_grid_spans", nothing)
   c === nothing || return c::NTuple{2,Matrix{Float32}}
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   hav(a1, o1, a2, o2) = 2 * 6.371e6 * asin(sqrt(sind((a2 - a1) / 2)^2 +
                                                cosd(a1) * cosd(a2) * sind((o2 - o1) / 2)^2))
   sx = zeros(Float32, Ni, Nj); sy = zeros(Float32, Ni, Nj)    # distance between the two stencil points
   for j in 1:Nj, i in 1:Ni
      i0 = max(i - 1, 1); i1 = min(i + 1, Ni); j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      sx[i, j] = hav(lat[i0, j], lon[i0, j], lat[i1, j], lon[i1, j])
      sy[i, j] = hav(lat[i, j0], lon[i, j0], lat[i, j1], lon[i, j1])
   end
   metconfig["_grid_spans"] = (sx, sy)
   return (sx, sy)
end

# On entry lev[:,:,k,indt] (k ≥ 2) holds the physical w [m/s] at slot k; on
# exit the level-index rate of the coordinate velocity (lev[:,:,1] untouched).
function _terrain_following_levdot!(metconfig, lev, zl, zs::AbstractMatrix, U, V, Ni, Nj, Nk, indt)
   sx, sy = _grid_spans_m!(metconfig, Ni, Nj)
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      @inbounds for i in 1:Ni
         i0 = max(i - 1, 1); i1 = min(i + 1, Ni)
         for k in 2:Nk
            zx = ((zl[i1,j,k,indt] + zs[i1,j]) - (zl[i0,j,k,indt] + zs[i0,j])) / max(sx[i,j], 1.0f0)
            zy = ((zl[i,j1,k,indt] + zs[i,j1]) - (zl[i,j0,k,indt] + zs[i,j0])) / max(sy[i,j], 1.0f0)
            lev[i,j,k,indt] -= U[i,j,k,indt] * zx + V[i,j,k,indt] * zy
         end
         for k in 2:Nk-1
            lev[i,j,k,indt] = levrate_from_w(lev[i,j,k,indt], zl[i,j,k-1,indt], zl[i,j,k+1,indt])
         end
         lev[i,j,Nk,indt] /= max(zl[i,j,Nk,indt] - zl[i,j,Nk-1,indt], 1.0f0)
      end
   end
   return nothing
end

# STACKED layout (2026-10-22).  The transport code reads zlev and levdot on HALF
# levels and u, v, T, p on the FULL levels stacked between them: interpolation4Dsingle
# reads a full-level field at index z + 0.5 (z ≥ 1.5; 2z − 1 below), the ERA5
# reader's layout (zlev = interface heights, slot s+1 = full level s).  On entry
# lev[:,:,s] holds the physical w on half level s, zl the half-level heights above
# ground, U/V the full-level winds at the mass points (slot s = the full level
# between zl[s-1] and zl[s]).  On exit the level-index rate of the coordinate
# velocity (w − u ∂z_s/∂x − v ∂z_s/∂y)/(∂z/∂k) on the half levels, with u, v on
# the half level = mean of the full levels above and below it; 0 at the ground
# (the half level of the terrain is a coordinate surface) and at the top.
function _terrain_following_levdot_stacked!(metconfig, lev, zl, zs::AbstractMatrix, U, V, Ni, Nj, Nk, indt)
   sx, sy = _grid_spans_m!(metconfig, Ni, Nj)
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      @inbounds for i in 1:Ni
         i0 = max(i - 1, 1); i1 = min(i + 1, Ni)
         lev[i,j,1,indt] = 0.0f0
         for s in 2:Nk-1
            zx = ((zl[i1,j,s,indt] + zs[i1,j]) - (zl[i0,j,s,indt] + zs[i0,j])) / max(sx[i,j], 1.0f0)
            zy = ((zl[i,j1,s,indt] + zs[i,j1]) - (zl[i,j0,s,indt] + zs[i,j0])) / max(sy[i,j], 1.0f0)
            uh = 0.5f0 * (U[i,j,s,indt] + U[i,j,s+1,indt])
            vh = 0.5f0 * (V[i,j,s,indt] + V[i,j,s+1,indt])
            lev[i,j,s,indt] = levrate_from_w(lev[i,j,s,indt] - uh * zx - vh * zy,
                                             zl[i,j,s-1,indt], zl[i,j,s+1,indt])
         end
         lev[i,j,Nk,indt] = 0.0f0
      end
   end
   return nothing
end

# One-cell distances on the sphere between mass points i−1 and i (dxw[i]) and
# j−1 and j (dyw[j]), cached per nest (C-grid stencil of the slope term below).
function _grid_spans1_m!(metconfig, Ni, Nj)
   c = get(metconfig, "_grid_spans1", nothing)
   c === nothing || return c::NTuple{2,Matrix{Float32}}
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   hav(a1, o1, a2, o2) = 2 * 6.371e6 * asin(sqrt(sind((a2 - a1) / 2)^2 +
                                                cosd(a1) * cosd(a2) * sind((o2 - o1) / 2)^2))
   dxw = zeros(Float32, Ni, Nj); dyw = zeros(Float32, Ni, Nj)
   for j in 1:Nj, i in 1:Ni
      i0 = max(i - 1, 1); j0 = max(j - 1, 1)
      dxw[i, j] = i0 == i ? hav(lat[i, j], lon[i, j], lat[i+1, j], lon[i+1, j]) : hav(lat[i0, j], lon[i0, j], lat[i, j], lon[i, j])
      dyw[i, j] = j0 == j ? hav(lat[i, j], lon[i, j], lat[i, j+1], lon[i, j+1]) : hav(lat[i, j0], lon[i, j0], lat[i, j], lon[i, j])
   end
   metconfig["_grid_spans1"] = (dxw, dyw)
   return (dxw, dyw)
end

# Coordinate velocity on terrain-following HALF levels from w, with the slope
# term in the C-grid form the driving model uses for its own kinematic w
# (2026-10-22): at a mass column, u·∂z/∂x = ½[u_w (z_i − z_{i−1})/Δx_w +
# u_e (z_{i+1} − z_i)/Δx_e] with the face winds u_w (west face, = Uf[i]) and u_e
# (east face, = Uf[i+1]) on the half level (mean of the two full levels around
# it), and likewise in y.  The centred form ½(u_w + u_e)(z_{i+1} − z_{i−1})/(2Δx)
# differs by ¼(u_w − u_e)(Δ₋ − Δ₊)/Δx, negative in valleys with converging
# drainage flow: there it gave a downward coordinate velocity at the first
# level and parcels were trapped at the ground (MesoNH, La Réunion: χ 2.5-2.9 in
# the lowest 100 m of the steepest blocks after 3 h without turbulence).
# Uf: (Ni or Ni+1, Nj, Nk) face winds on the full-level slots, face i = west face
# of mass column i (the last column reuses its west face when there is no Ni+1);
# Vf likewise on the south faces.  lev on entry: w on the half levels.
function _terrain_following_levdot_stacked_cgrid!(metconfig, lev, zl, zs::AbstractMatrix, Uf, Vf, Ni, Nj, Nk, indt)
   dxw, dyw = _grid_spans1_m!(metconfig, Ni, Nj)
   nfu = size(Uf, 1); nfv = size(Vf, 2)
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         i0 = max(i - 1, 1); i1 = min(i + 1, Ni); j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
         ie = min(i + 1, nfu); jn = min(j + 1, nfv)
         dxe = i1 == i ? dxw[i, j] : dxw[i1, j]; dyn = j1 == j ? dyw[i, j] : dyw[i, j1]
         lev[i,j,1,indt] = 0.0f0
         for s in 2:Nk-1
            zc = zl[i,j,s,indt] + zs[i,j]
            uw = 0.5f0 * (Uf[i, j,s,indt] + Uf[i, j,s+1,indt]); ue = 0.5f0 * (Uf[ie,j,s,indt] + Uf[ie,j,s+1,indt])
            vs = 0.5f0 * (Vf[i,j, s,indt] + Vf[i,j, s+1,indt]); vn = 0.5f0 * (Vf[i,jn,s,indt] + Vf[i,jn,s+1,indt])
            cx = 0.5f0 * (uw * (zc - (zl[i0,j,s,indt] + zs[i0,j])) / max(dxw[i,j], 1.0f0) +
                          ue * ((zl[i1,j,s,indt] + zs[i1,j]) - zc) / max(dxe, 1.0f0))
            cy = 0.5f0 * (vs * (zc - (zl[i,j0,s,indt] + zs[i,j0])) / max(dyw[i,j], 1.0f0) +
                          vn * ((zl[i,j1,s,indt] + zs[i,j1]) - zc) / max(dyn, 1.0f0))
            lev[i,j,s,indt] = levrate_from_w(lev[i,j,s,indt] - cx - cy, zl[i,j,s-1,indt], zl[i,j,s+1,indt])
         end
         lev[i,j,Nk,indt] = 0.0f0
      end
   end
   return nothing
end

# ── WRF (NetCDF wrfout, Arakawa-C grid, terrain-following eta) ───────────────
# Rewritten 2026-10-13 (reported by Yue Jia; see CHANGES.md).  The previous
# path (checkmetfields!) wrote only rho back: T stayed WRF's perturbation
# potential temperature θ' = θ − 300 K, theta stayed zero, puv was never
# filled; heights were the eta INTERFACES while T, p, u, v are on the mass
# levels between them (half a layer apart, the top mass level left empty);
# U and V were used on their staggered faces as if they were mass points; and
# no surface-turbulence derivation ran, so tke stayed zero (σw floors).
#
# Loader conventions on entry (_read_netcdf_file!, metconfig Nk = interfaces):
#   mass fields (T, QVAPOR, P, PB)  at k = 1 … Nk−1  (mass level k)
#   interface fields (PH, PHB, W)   at k = 1 … Nk    (interface k = bottom of mass level k)
#   U on Ni+1 x-faces, V on Nj+1 y-faces, GRID-relative (what the Lambert
#   advection needs: displacements are along the grid axes, with the map factor)
# Output: the common layout — k = 1 surface slot (zlev = 0, p = PSFC, T = T2
# when mapped), k = 2 … Nk the mass levels at their own heights above ground,
# u/v destaggered to the mass points, levdot from W destaggered to the mass
# levels, rho/theta/drhodz from p, T, q, and the Hanna tke from UST, PBLH and
# HFX (positive upward in WRF, which is the convention _derive_surface_params!
# expects).
const _WRF_RCP = 2.0f0 / 7.0f0          # WRF's R_d/c_p (r_d = 287, cp = 7 r_d / 2)

function metwrffieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   avg = get(metconfig, "wrf_winds", "instant") == "avgflx"
   any(x -> occursin("AVGFLX", string(x)), metconfig["input_vars_3D"]) && !avg &&
      error("WRF: AVGFLX fields are mapped but the namelist does not declare them — use " *
            "listvariables/WRFv2_avgflx.vars (#info wrf_winds = avgflx)")
   if avg
      for key in ("mu", "mub")
         haskey(metfields2D, key) ||
            error("WRF avgflx: '$key' missing — WRFv2_avgflx.vars maps mu = MU, mub = MUB")
      end
      for key in ("wrf_mapfac_uy", "wrf_mapfac_vx", "wrf_mapfac_my", "wrf_znu")
         haskey(metconfig, key) ||
            error("WRF avgflx: $(uppercase(key[5:end])) not found in the wrfout file")
      end
   end
   for key in ("Ux", "Vy", "T", "p", "pb", "ph", "phb")
      haskey(metfields3D, key) ||
         error("WRF reader: '$key' missing — WRFv2.vars maps Ux=U, Vy=V, W=W, T=T, " *
               "p=P, pb=PB, ph=PH, phb=PHB")
   end
   Ni = Int(metconfig["Ni"]); Nj = Int(metconfig["Nj"]); Nk = Int(metconfig["Nk"])
   Km = Nk - 1                                            # mass levels
   T  = metfields3D["T"];  zl = metfields3D["zlev"]; pu = metfields3D["puv"]
   pp = metfields3D["p"];  pb = metfields3D["pb"]
   ph = metfields3D["ph"]; phb = metfields3D["phb"]
   U  = metfields3D["Ux"]; V = metfields3D["Vy"]
   qv = haskey(metfields3D, "qv") ? metfields3D["qv"] : nothing
   W  = haskey(metfields3D, "W") ? metfields3D["W"] : nothing
   lev = metfields3D["levdot"]
   # positions read from the file at setup (FieldPos, initialize.jl): U, V on the
   # faces or at the mass points, W and tke on the w-levels or the mass levels
   fp = get(metconfig, "field_pos", Dict{String,FieldPos}())
   uedge = haskey(fp, "Ux") ? fp["Ux"].x == :face : size(U, 1) == Ni + 1
   vedge = haskey(fp, "Vy") ? fp["Vy"].y == :face : size(V, 2) == Nj + 1
   (uedge == (size(U, 1) == Ni + 1) && vedge == (size(V, 2) == Nj + 1)) ||
      error("WRF: the U/V arrays ($(size(U, 1))×$(size(V, 2))) do not match the positions read from the file (x-faces $uedge, y-faces $vedge)")
   wfull = W !== nothing && haskey(fp, "W") && fp["W"].z == :full
   avg && wfull && error("WRF avgflx: AVGFLX_WWM must be on the w-levels (bottom_top_stag)")
   has2t = haskey(metfields2D, "2t"); hassp = haskey(metfields2D, "sp")
   # raw staggered winds of this slot: rows/columns are rewritten in place below
   # (by other threads too), so destagger from a snapshot
   Ur = U[:, :, 1:Km, indt]; Vr = V[:, :, 1:Km, indt]
   # Time-averaged mass-coupled fluxes (2026-10-21).  AVGFLX_RUM = (c1h μ_u + c2h)
   # u / m_y, AVGFLX_RVM = (c1h μ_v + c2h) v / m_x, AVGFLX_WWM = Ω = μ η̇ / m_y
   # (μ = MU + MUB the dry column mass; verified on a WRF-Chem 4.2.2 file: the
   # column integral of the flux divergence reproduces WWM with slope 1.000 and
   # closes at the top), each the mean over the history interval ending at the
   # output time.  Being WRF's own mass fluxes they are continuity-consistent,
   # which the instantaneous U, V, W snapshots are not (the coordinate velocity
   # derived from W was 3× noisier than η̇ in the lowest 50 m).
   if avg
      (uedge && vedge) || error("WRF avgflx: AVGFLX_RUM/RVM must be on the staggered faces")
      # WRF writes the first average one history interval after its start: the
      # AVGFLX fields of the start time are identically zero (measured on the
      # WRF-Chem sample: max |RUM| 0 at 00:00, 1.8e6 at 01:00) — no silent calm
      maximum(abs, view(U, :, :, 1:Km, indt)) == 0 &&
         error("WRF avgflx: AVGFLX_RUM is zero at this output time (the WRF simulation start, " *
               "before the first averaging interval) — start the EASLEE run one history " *
               "interval later, or read the instantaneous winds with WRFv2.vars")
      mud  = metfields2D["mu"][1:Ni, 1:Nj, indt] .+ metfields2D["mub"][1:Ni, 1:Nj, indt]
      c1h  = metconfig["wrf_c1h"]::Vector{Float32}; c2h = metconfig["wrf_c2h"]::Vector{Float32}
      c1f  = metconfig["wrf_c1f"]::Vector{Float32}; c2f = metconfig["wrf_c2f"]::Vector{Float32}
      mfuy = metconfig["wrf_mapfac_uy"]::Matrix{Float32}; mfvx = metconfig["wrf_mapfac_vx"]::Matrix{Float32}
      mfmy = metconfig["wrf_mapfac_my"]::Matrix{Float32}
      Threads.@threads :static for m in 1:Km
         @inbounds for j in 1:Nj, i in 1:Ni+1
            mu_u = 0.5f0 * (mud[clamp(i - 1, 1, Ni), j] + mud[min(i, Ni), j])
            Ur[i,j,m] *= mfuy[i,j] / (c1h[m] * mu_u + c2h[m])
         end
         @inbounds for j in 1:Nj+1, i in 1:Ni
            mu_v = 0.5f0 * (mud[i, clamp(j - 1, 1, Nj)] + mud[i, min(j, Nj)])
            Vr[i,j,m] *= mfvx[i,j] / (c1h[m] * mu_v + c2h[m])
         end
      end
      # dη/dk on the half levels (slot s = w-level s, η = ZNW[s])
      znw = metconfig["wrf_znw"]::Vector{Float32}
      detadk = [s == 1 || s == Nk ? 1.0f0 : 0.5f0 * (znw[s+1] - znw[s-1]) for s in 1:Nk]
   end
   zsfc = zeros(Float32, Ni, Nj)                          # terrain height ASL (PH+PHB at the ground)
   Threads.@threads :static for j in 1:Nj
      zw = Vector{Float32}(undef, Nk)
      @inbounds for i in 1:Ni
         for k in 1:Nk                                    # interface heights, 1 = the ground
            zw[k] = (ph[i,j,k,indt] + phb[i,j,k,indt]) / _G0
         end
         zs = zw[1]; zsfc[i,j] = zs
         # mass level m → slot m + 1, top down: slot m still holds level m when read
         for m in Km:-1:1
            k = m + 1
            P  = pp[i,j,m,indt] + pb[i,j,m,indt]
            T[i,j,k,indt]  = (T[i,j,m,indt] + 300.0f0) * (P / 100000.0f0)^_WRF_RCP
            pu[i,j,k,indt] = P
            if qv !== nothing
               r = max(qv[i,j,m,indt], 0.0f0)
               qv[i,j,k,indt] = r / (1.0f0 + r)            # QVAPOR is a mixing ratio
            end
            U[i,j,k,indt] = uedge ? 0.5f0 * (Ur[i,j,m] + Ur[i+1,j,m]) : Ur[i,j,m]
            V[i,j,k,indt] = vedge ? 0.5f0 * (Vr[i,j,m] + Vr[i,j+1,m]) : Vr[i,j,m]
         end
         # half levels = the w-levels: zlev and the vertical motion (stacked layout)
         for s in 1:Nk
            zl[i,j,s,indt] = zw[s] - zs
            if W !== nothing
               lev[i,j,s,indt] = avg ?                     # η̇ = Ω m_y / (c1f μ + c2f)
                  (s == 1 || s == Nk ? 0.0f0 : W[i,j,s,indt] * mfmy[i,j] / (c1f[s] * mud[i,j] + c2f[s]) / detadk[s]) :
                  !wfull ? W[i,j,s,indt] :                 # physical w on the w-levels, converted below
                  (s == 1 || s == Nk ? 0.0f0 :             # w on the mass levels (slots 1…Km as loaded):
                   0.5f0 * (W[i,j,s-1,indt] + W[i,j,s,indt]))   # w-level s between mass levels s−1 and s
            end
         end
         # surface slot
         zl[i,j,1,indt] = 0.0f0
         pu[i,j,1,indt] = hassp ? metfields2D["sp"][i,j,indt] :          # full level 1 at ~zl[2]/2
                          pu[i,j,2,indt] * exp(_G0 * 0.5f0 * zl[i,j,2,indt] / (_RD * T[i,j,2,indt]))
         T[i,j,1,indt]  = has2t ? metfields2D["2t"][i,j,indt] : T[i,j,2,indt]
         qv === nothing || (qv[i,j,1,indt] = qv[i,j,2,indt])
         # surface slot wind = the lowest level's, as in every other reader: below
         # the lowest level the Monin-Obukhov option ([physics] surface_layer_wind,
         # on by default) rescales it onto the similarity profile.  U10/V10 here
         # (c6386f7) made that reduction apply on top of a U10 blend (2026-10-21).
         U[i,j,1,indt]  = U[i,j,2,indt]
         V[i,j,1,indt]  = V[i,j,2,indt]
         lev[i,j,1,indt] = 0.0f0
      end
   end
   # the extra face column/row of the staggered arrays: copy the last mass point
   uedge && (@views U[Ni+1, :, :, indt] .= U[Ni, :, :, indt])
   vedge && (@views V[:, Nj+1, :, indt] .= V[:, Nj, :, indt])
   # levdot: w [m/s] on the w-levels (= the half levels) → coordinate velocity → level-index rate
   if W === nothing
      fill!(view(lev, :, :, :, indt), 0.0f0)
   elseif !avg                    # (avgflx: η̇ is already the coordinate velocity)
      # centred slope term: measured better than the C-grid form on the WRF-Chem
      # sample (0-100 m χ 1.043 vs 1.094 after 4 h; the C-grid form is used for MesoNH)
      _terrain_following_levdot_stacked!(metconfig, lev, zl, zsfc, U, V, Ni, Nj, Nk, indt)
   else
      # WRF never computes the averaged Ω on the outermost row and column of mass
      # points (measured: column continuity off by 19-75 Pa/s there, exact to
      # 1e-6 everywhere else) — take the adjacent interior cell's rate
      @views begin
         lev[1,  :, 2:Nk, indt] .= lev[2,    :, 2:Nk, indt]
         lev[Ni, :, 2:Nk, indt] .= lev[Ni-1, :, 2:Nk, indt]
         lev[:,  1, 2:Nk, indt] .= lev[:,    2, 2:Nk, indt]
         lev[:, Nj, 2:Nk, indt] .= lev[:, Nj-1, 2:Nk, indt]
      end
   end

   _thermo_from_pTq!(metfields3D, Nk, indt)
   haskey(fp, "tke") && @warn "WRF: '$(fp["tke"].input)' is mapped to tke but the WRF reader replaces tke with the Hanna fill (from UST, PBLH, HFX); the archived TKE is not used" maxlog=1
   _surface_turbulence!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   return nothing
end

# OMEGA = Dp/Dt on hybrid sigma-pressure levels → level-index rate.  A parcel
# keeps its level index while it moves horizontally, so along the level its
# pressure already changes by u·∂p_k/∂x + v·∂p_k/∂y; the level rate must come
# from the remainder ω − u·∂p_k/∂x − v·∂p_k/∂y (the tendency ∂p_k/∂t =
# b_k ∂p_s/∂t is neglected).  Until 2026-10-13 ω was used whole: measured on
# MERRA-2 2022-01-18 over terrain above 1000 m, the lowest level's ω (0.198 Pa/s
# mean |ω|) correlates at 0.95 with V·∇p_k and the coordinate part is 0.067 —
# the vertical motion near the ground over mountains was ~3× too large.
# Regular lon/lat grid (MERRA-2), periodic in longitude when global.
function _hybrid_p_levdot!(metconfig, lev, p, U, V, Nj, Nk, indt, globalx::Bool)
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   N0 = length(lon); R = 6.371e6; dl = deg2rad(lon[2] - lon[1])
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      cl = cosd(lat[j]); spy = Float32(R * deg2rad(lat[j1] - lat[j0]))
      @inbounds for i in 1:N0
         if globalx
            i0 = i == 1 ? N0 : i - 1; i1 = i == N0 ? 1 : i + 1; spx = Float32(2 * dl * R * cl)
         else
            i0 = max(i - 1, 1); i1 = min(i + 1, N0); spx = Float32((i1 - i0) * dl * R * cl)
         end
         for k in 2:Nk
            px = cl > 1.0e-3 ? (p[i1,j,k,indt] - p[i0,j,k,indt]) / spx : 0.0f0
            py = (p[i,j1,k,indt] - p[i,j0,k,indt]) / spy
            lev[i,j,k,indt] -= U[i,j,k,indt] * px + V[i,j,k,indt] * py
         end
         for k in 2:Nk-1
            lev[i,j,k,indt] *= 2.0f0 / min(p[i,j,k+1,indt] - p[i,j,k-1,indt], -1.0f0)
         end
      end
      if globalx && size(lev, 1) > N0                  # wrap column
         @views lev[N0+1,j,:,indt] .= lev[1,j,:,indt]
      end
   end
   return nothing
end

# STACKED layout helpers for readers whose heights and vertical velocity come on
# the FULL levels (MERRA-2 H/OMEGA, GFS isobaric gh/w, AROME-HDF fixed heights):
# the transport reads zlev and levdot on the half levels between them (see
# _terrain_following_levdot_stacked!).
# Full-level heights above ground in zl (slot 1 = the ground) → half-level
# heights in place: zlev[s] = midpoint of the full levels s and s+1, 0 at the
# ground, the top half a layer above the top full level.
function _stack_half_levels!(zl, Ni, Nj, Nk, indt)
   Threads.@threads :static for j in 1:Nj
      @inbounds for i in 1:Ni
         ztop = zl[i,j,Nk,indt] + 0.5f0 * (zl[i,j,Nk,indt] - zl[i,j,Nk-1,indt])
         zl[i,j,1,indt] = 0.0f0
         for s in 2:Nk-1
            zl[i,j,s,indt] = 0.5f0 * (zl[i,j,s,indt] + zl[i,j,s+1,indt])
         end
         zl[i,j,Nk,indt] = ztop
      end
   end
   return nothing
end
# ω [Pa/s] on the full levels (lev[:,:,s], p the full-level pressures) → level-
# index rate of the coordinate velocity on the half levels: on half level s
# (between the full levels s and s+1) ω_h − u_h ∂p_h/∂x − v_h ∂p_h/∂y over
# ∂p/∂k = p[s+1] − p[s] (one index between the two full levels), ω_h, u_h, v_h,
# p_h the means of the two full levels.  0 at the ground (slot 1) and the top.
# Ascending s in place: lev[s+1] still holds a full-level ω when read.
function _hybrid_p_levdot_stacked!(metconfig, lev, p, U, V, Nj, Nk, indt, globalx::Bool)
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   ndims(lat) == 2 && (lat = lat[1, :]); ndims(lon) == 2 && (lon = lon[:, 1])   # GRIB: 2-D
   N0 = length(lon); R = 6.371e6; dl = deg2rad(lon[2] - lon[1])
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      cl = cosd(lat[j]); spy = Float32(R * deg2rad(lat[j1] - lat[j0]))
      @inbounds for i in 1:N0
         if globalx
            i0 = i == 1 ? N0 : i - 1; i1 = i == N0 ? 1 : i + 1; spx = Float32(2 * dl * R * cl)
         else
            i0 = max(i - 1, 1); i1 = min(i + 1, N0); spx = Float32((i1 - i0) * dl * R * cl)
         end
         lev[i,j,1,indt] = 0.0f0
         for s in 2:Nk-1
            px = cl > 1.0e-3 ? 0.5f0 * ((p[i1,j,s,indt] + p[i1,j,s+1,indt]) - (p[i0,j,s,indt] + p[i0,j,s+1,indt])) / spx : 0.0f0
            py = 0.5f0 * ((p[i,j1,s,indt] + p[i,j1,s+1,indt]) - (p[i,j0,s,indt] + p[i,j0,s+1,indt])) / spy
            wh = 0.5f0 * (lev[i,j,s,indt] + lev[i,j,s+1,indt])
            uh = 0.5f0 * (U[i,j,s,indt] + U[i,j,s+1,indt]); vh = 0.5f0 * (V[i,j,s,indt] + V[i,j,s+1,indt])
            lev[i,j,s,indt] = (wh - uh * px - vh * py) / min(p[i,j,s+1,indt] - p[i,j,s,indt], -1.0f0)
         end
         lev[i,j,Nk,indt] = 0.0f0
      end
      if globalx && size(lev, 1) > N0                  # wrap column
         @views lev[N0+1,j,:,indt] .= lev[1,j,:,indt]
      end
   end
   return nothing
end

# ── MERRA2 (NetCDF, GEOS eta levels; PL and H archived directly) ─────────────
function metmerra2fieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   Ni = Int(metconfig["Ni"]) + (globalx ? 1 : 0)
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   _shift_up_one!(metfields3D, ("Ux","Vy","T","qv","levdot","puv","zlev","cc","kh"), Nk, indt)

   sp   = metfields2D["sp"][:,:,indt]
   zter = metfields2D["phis"][:,:,indt] ./ _G0      # surface elevation [m]
   zl = metfields3D["zlev"]
   for k in 2:Nk
      zl[:,:,k,indt] = max.(zl[:,:,k,indt] .- zter, Float32(k - 1))
   end
   zl[:,:,1,indt] .= 0.0f0

   metfields3D["puv"][:,:,1,indt] = sp
   for key in ("Ux","Vy","T","qv","cc")
      haskey(metfields3D, key) || continue
      metfields3D[key][:,:,1,indt] = metfields3D[key][:,:,2,indt]
   end

   # STACKED layout (2026-10-22): H is archived on the mid-layer (full) levels
   # with U, V, T, PL — zlev goes to the half levels between them, and OMEGA to
   # the coordinate velocity on those half levels (_hybrid_p_levdot_stacked!).
   _stack_half_levels!(zl, Ni, Nj, Nk, indt)
   _thermo_from_pTq!(metfields3D, Nk, indt)
   _hybrid_p_levdot_stacked!(metconfig, metfields3D["levdot"], metfields3D["puv"],
                             metfields3D["Ux"], metfields3D["Vy"], Nj, Nk, indt, globalx)

   # Edge (layer-interface) fields → mid-layer slots.  The mst_Ne collection is
   # on the 73 layer EDGES, loaded bottom-up with edge 1 = surface, so slot k
   # held the TOP edge of the layer slot k represents.  Everything downstream
   # treats umf and tpf as mid-layer values: the convection scheme builds its
   # interface flux as 0.5·(umf[k] + umf[k+1]) and wet deposition interpolates
   # tpf at the parcel's level like T or q.  Used raw, both sat half a layer
   # too high (2026-10-13).  Slot 1 keeps the surface edge (umf = 0 there, tpf
   # = surface precipitation); slot k ≥ 2 = mean of its bottom and top edges.
   for key in ("umf","pflcu","pficu","pflls","pfils")
      haskey(metfields3D, key) || continue
      A = metfields3D[key]
      for k in Nk:-1:2                       # descending: A[k-1] still holds an edge
         @views A[:,:,k,indt] .= 0.5f0 .* (A[:,:,k-1,indt] .+ A[:,:,k,indt])
      end
   end

   # Total precipitation flux for wet deposition: sum of the mst_Ne fluxes
   # (convective + large-scale/anvil, rain + ice), now on mid-layer slots.
   if haskey(metfields3D, "tpf")
      tpf = metfields3D["tpf"]
      fill!(view(tpf, :, :, :, indt), 0.0f0)
      for key in ("pflcu","pficu","pflls","pfils")
         haskey(metfields3D, key) || continue
         tpf[:,:,:,indt] .+= max.(metfields3D[key][:,:,:,indt], 0.0f0)
      end
   end

   _surface_turbulence!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   return nothing
end

# Kinematic ω at the ground, V·∇p_s, with the wind of the first isobaric level
# above ground (∂p_s/∂t is not available from one file and is neglected).
function _gfs_surface_omega(metconfig, sp, metfields3D, Ni, Nj, plev, indt)
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   ndims(lat) == 2 && (lat = lat[1, :]); ndims(lon) == 2 && (lon = lon[:, 1])
   N0 = length(lon); R = 6.371e6; dl = deg2rad(lon[2] - lon[1])
   U = metfields3D["Ux"]; V = metfields3D["Vy"]; Nk = length(plev) + 1
   oms = zeros(Float32, Ni, Nj)
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      cl = cosd(lat[j]); spy = Float32(R * deg2rad(lat[j1] - lat[j0]))
      @inbounds for i in 1:min(Ni, N0)
         i0 = i == 1 ? (Ni > N0 ? N0 : 1) : i - 1
         i1 = i == N0 ? (Ni > N0 ? 1 : N0) : i + 1
         spx = Float32(R * cl * dl * (Ni > N0 ? 2 : (i1 - i0)))
         kg = 2
         while kg < Nk && plev[kg-1] >= sp[i,j]
            kg += 1
         end
         px = cl > 1.0e-3 ? (sp[i1,j] - sp[i0,j]) / spx : 0.0f0
         py = (sp[i,j1] - sp[i,j0]) / spy
         oms[i,j] = U[i,j,kg,indt] * px + V[i,j,kg,indt] * py
      end
      Ni > N0 && (oms[Ni,j] = oms[1,j])
   end
   return oms
end

# GFS isobaric levels → a hybrid σ–p coordinate (2026-10-22).  The isobaric
# ladder runs below the terrain; the 2026-10-13 surface ramp gave levels whose
# definition jumped between columns wherever a pressure level crossed the
# ground, and the transport piled air up at the ground there (χ 1.29 in the
# lowest 100 m after 4 h on GFS 0.25°, steep blocks 1.34).  Now slot s ≥ 2 sits
# at p_s = p_ref,s − B_s·(p0 − p_sfc), p_ref the isobaric ladder, B linear in
# p_ref from 0 at 200 hPa (pure pressure above) to 1 at p0: terrain-following at
# the ground, GFS's own levels aloft, and monotone for any p_sfc > 200 hPa
# (∂p_s/∂p_ref = (p_sfc − 200 hPa)/(p0 − 200 hPa) > 0).  The first version
# (p_ref·(1 − w + w·p_sfc/p0), w linear from 1 at 850 to 0 at 200 hPa) folded
# where p_sfc < ~560 hPa: over Tibet p(800) > p(850), two half levels 1 m apart
# and a level rate of 2.6 index/s — parcels jumped 17 levels in one sub-step.  u, v, T, q, cc, ω and the height above ground are
# interpolated in ln p from the levels above the ground and the surface values
# (T = 2 m T, u, v, q, cc of the first level above ground, ω = V·∇p_s, z = 0).
function _gfs_to_hybrid!(metfields3D, metfields2D, plev, Ni, Nj, Nk, indt, oms::AbstractMatrix)
   p0 = 101325.0f0; p_prs = 20000.0f0
   bv = [clamp((plev[s-1] - p_prs) / (p0 - p_prs), 0.0f0, 1.0f0) for s in 2:Nk]
   names = [k for k in ("Ux","Vy","T","qv","cc","levdot","zlev") if haskey(metfields3D, k)]
   arrs  = [metfields3D[k] for k in names]
   p = metfields3D["puv"]; T2 = haskey(metfields2D, "2t") ? metfields2D["2t"] : nothing
   iT = findfirst(==("T"), names); iw = findfirst(==("levdot"), names); iz = findfirst(==("zlev"), names)
   Threads.@threads :static for j in 1:Nj
      nf = length(arrs)
      lnps = Vector{Float32}(undef, Nk); src = Matrix{Float32}(undef, nf, Nk); out = Matrix{Float32}(undef, nf, Nk)
      ptg = Vector{Float32}(undef, Nk)
      @inbounds for i in 1:Ni
         sp = p[i,j,1,indt]
         kg = 2
         while kg < Nk && plev[kg-1] >= sp
            kg += 1
         end
         # source column: the surface, then the isobaric levels above it
         n = 1; lnps[1] = log(sp)
         for (a, A) in enumerate(arrs)
            src[a, 1] = A[i,j,kg,indt]                       # u, v, q, cc of the first level above ground
         end
         src[iT, 1] = T2 === nothing ? arrs[iT][i,j,kg,indt] : T2[i,j,indt]
         src[iw, 1] = oms[i,j]
         src[iz, 1] = 0.0f0
         for k in kg:Nk
            n += 1; lnps[n] = log(plev[k-1])
            for (a, A) in enumerate(arrs)
               src[a, n] = A[i,j,k,indt]
            end
         end
         # target pressures and ln-p interpolation (lnps descending)
         ptg[1] = sp
         for s in 2:Nk
            ptg[s] = plev[s-1] - bv[s-1] * (p0 - sp)
         end
         for s in 1:Nk
            lp = log(ptg[s]); m = 1
            while m < n - 1 && lnps[m+1] > lp
               m += 1
            end
            # a surface pressure within Float32 resolution of an isobaric level
            # (97500.055 Pa over 975 hPa) gives a zero-width interval: 0/0 made
            # every full-level surface value NaN in 80 columns of a GFS file
            d = lnps[m+1] - lnps[m]
            f = d != 0.0f0 ? clamp((lp - lnps[m]) / d, 0.0f0, 1.0f0) : 0.0f0
            for a in 1:nf
               out[a, s] = src[a, m] + f * (src[a, m+1] - src[a, m])
            end
         end
         for s in 1:Nk
            p[i,j,s,indt] = ptg[s]
            for (a, A) in enumerate(arrs)
               A[i,j,s,indt] = out[a, s]
            end
         end
         arrs[iz][i,j,1,indt] = 0.0f0
         for s in 2:Nk                                       # keep the heights strictly increasing
            arrs[iz][i,j,s,indt] = max(arrs[iz][i,j,s,indt], arrs[iz][i,j,s-1,indt] + 1.0f0)
         end
      end
   end
   return nothing
end

# ── NCEP GFS (GRIB2, isobaric levels) ─────────────────────────────────────────
function metgfsfieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   get(metconfig, "_layout_logged", false) || (metconfig["_layout_logged"] = true;
      @info "met field positions (GFS isobaric GRIB): every field on the regular grid at the same points (no staggering); " *
            "the isobaric levels are remapped to hybrid σ–p full levels (slots 2…Nk, slot 1 = surface values); " *
            "zlev = midpoints of the full-level heights (half levels); ω → coordinate velocity on the half levels")
   Ni = Int(metconfig["Ni"]) + (globalx ? 1 : 0)
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   plev = Float32.(metconfig["levelvalues"]) .* 100.0f0   # hPa → Pa, descending

   _shift_up_one!(metfields3D, ("Ux","Vy","T","qv","levdot","zlev","cc"), Nk, indt)

   sp   = metfields2D["sp"][:,:,indt]
   orog = haskey(metfields2D, "orog") ? metfields2D["orog"][:,:,indt] : zeros(Float32, Ni, Nj)

   # Heights above ground from geopotential height (levels below the terrain
   # are not used: _gfs_to_hybrid! interpolates from the levels above it)
   zl = metfields3D["zlev"]
   zl[:,:,1,indt] .= 0.0f0
   for k in 2:Nk
      zl[:,:,k,indt] = zl[:,:,k,indt] .- orog
   end

   p = metfields3D["puv"]
   p[:,:,1,indt] = sp
   for k in 2:Nk
      p[:,:,k,indt] .= plev[k-1]
   end

   if haskey(metfields2D, "2t")
      metfields3D["T"][:,:,1,indt] = metfields2D["2t"][:,:,indt]
   else
      metfields3D["T"][:,:,1,indt] = metfields3D["T"][:,:,2,indt]
   end
   for key in ("Ux","Vy","qv","cc")
      haskey(metfields3D, key) || continue
      metfields3D[key][:,:,1,indt] = metfields3D[key][:,:,2,indt]
   end

   # Hybrid σ–p full levels (_gfs_to_hybrid!), then the STACKED layout
   # (2026-10-22): zlev on the half levels between the full levels and ω to the
   # coordinate velocity on those half levels (_hybrid_p_levdot_stacked!).  Until
   # 2026-10-22: isobaric slots, a surface ramp below the first level above
   # ground, levdot = 2ω/Δp with a fixed ladder spacing, zero across the ramp.
   oms = _gfs_surface_omega(metconfig, sp, metfields3D, Ni, Nj, plev, indt)
   _gfs_to_hybrid!(metfields3D, metfields2D, plev, Ni, Nj, Nk, indt, oms)
   _stack_half_levels!(zl, Ni, Nj, Nk, indt)
   _thermo_from_pTq!(metfields3D, Nk, indt)
   _hybrid_p_levdot_stacked!(metconfig, metfields3D["levdot"], p, metfields3D["Ux"], metfields3D["Vy"],
                             Nj, Nk, indt, globalx)

   # GFS cloud fraction is archived in percent → normalise to 0–1 (the wet
   # deposition thresholds and the cloud-top search expect fractions).
   if haskey(metfields3D, "cc")
      metfields3D["cc"][:,:,:,indt] = clamp.(metfields3D["cc"][:,:,:,indt] ./ 100.0f0, 0.0f0, 1.0f0)
   end
   # 3-D precipitation flux for wet deposition from the surface rate.
   _tpf_from_surface_precip!(metfields3D, metfields2D, indt, Ni, Nj, Nk)

   _surface_turbulence!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   return nothing
end

# Slots nhk+1 … k0−1 (isobaric levels not above the top height level nhk + 50 m)
# → heights evenly spaced between z(nhk) and z(k0), pressure log-linear, T/u/v/rh
# linear in height through the original profile (the block top, then the
# isobaric levels above it); k0 = the first isobaric level above that leaves at
# least _STITCH_DZMIN per ramped layer.  Until 2026-10-22 k0 was the first level
# above z(nhk) + 50 m and the ramped levels could be crushed into a few tens of
# metres: a 17 m layer at 3.8 km gave the continuity w (6.9 m/s there) a level
# rate of 0.4 index/s, i.e. an index-space step of ~24 levels in 60 s.  NaN
# columns (bitmap-missing) are left as they are.
const _STITCH_DZMIN = 100.0f0
function _arome_hdf_stitch_ramp!(metfields3D, nhk, Nk, Ni, Nj, indt)
   zl = metfields3D["zlev"]; p = metfields3D["puv"]
   arrs = [metfields3D[k] for k in ("T", "Ux", "Vy", "rh") if haskey(metfields3D, k)]
   na = length(arrs)
   Threads.@threads :static for j in 1:Nj
      zo = Vector{Float32}(undef, Nk); po = Vector{Float32}(undef, Nk); vo = Matrix{Float32}(undef, na, Nk)
      @inbounds for i in 1:Ni
         ztop = zl[i,j,nhk,indt]
         isfinite(ztop) || continue
         k0 = nhk + 1
         while k0 <= Nk && !(zl[i,j,k0,indt] > ztop + 50.0f0)      # NaN fails the test
            k0 += 1
         end
         k0 > Nk && continue
         while k0 < Nk && isfinite(zl[i,j,k0+1,indt]) &&
               zl[i,j,k0,indt] - ztop < _STITCH_DZMIN * Float32(k0 - nhk)
            k0 += 1
         end
         k0 == nhk + 1 && continue
         # original profile on [z(nhk), z(k0)], heights increasing
         n = 1; zo[1] = ztop; po[1] = p[i,j,nhk,indt]
         for a in 1:na; vo[a, 1] = arrs[a][i,j,nhk,indt]; end
         for k in nhk+1:k0
            z = zl[i,j,k,indt]
            (isfinite(z) && z > zo[n] + 50.0f0) || continue
            n += 1; zo[n] = z; po[n] = p[i,j,k,indt]
            for a in 1:na; vo[a, n] = arrs[a][i,j,k,indt]; end
         end
         z1 = zo[n]; nr = Float32(k0 - nhk)
         for k in nhk+1:k0-1
            zt = ztop + Float32(k - nhk) / nr * (z1 - ztop)
            m = 1
            while m < n - 1 && zo[m+1] < zt
               m += 1
            end
            f = clamp((zt - zo[m]) / max(zo[m+1] - zo[m], 1.0f-3), 0.0f0, 1.0f0)
            zl[i,j,k,indt] = zt
            p[i,j,k,indt]  = po[m] * (po[m+1] / po[m])^f
            for a in 1:na
               arrs[a][i,j,k,indt] = vo[a, m] + f * (vo[a, m+1] - vo[a, m])
            end
         end
      end
   end
   return nothing
end

# Vertical velocity diagnosed from mass continuity, for datasets that archive
# none (AROME-HDF), STACKED layout, regular lon/lat grid (2026-10-22).  Layer of
# full level s (s ≥ 2) between zlev[s-1] and zlev[s]: mass per area m_s = ρ_s Δz_s,
# horizontal divergence D_s of m_s·V_s (centred, on the sphere); the upward mass
# flux through half level s is
#     F_s = (M_{≤s}/M)·ΣD − Σ_{j≤s} D_j      (F = 0 at the ground and the top)
# with M_{≤s} the air mass below half level s and M the column's: the column's
# net divergence (its mass tendency) is shared out in proportion to mass, as in
# σ / hybrid coordinates (B·∂p_s/∂t).  The level rate is F_s / (ρ_h·(zf_{s+1} −
# zf_s)), ρ_h and zf the half-level density and the full-level heights.  Without
# it the dataset had levdot = 0 and the well-mixed test failed by construction
# (block spread 2-3× Poisson after 1 h).  The first version (680558d) integrated
# F_s = −Σ_{j≤s} D_j upward with nothing imposed at the top: the residual grew
# with height (p99.99 of |w| 5 m/s at 3-7 km, 11 m/s at 15 km, max 20 m/s) and
# drove parcels into the model top.
function _continuity_levdot!(metconfig, lev, zl, rho, U, V, Ni, Nj, Nk, indt)
   lat = Float64.(metconfig["latitudes"]); lon = Float64.(metconfig["longitudes"])
   ndims(lat) == 2 && (lat = lat[1, :]); ndims(lon) == 2 && (lon = lon[:, 1])
   R = 6.371e6; dl = deg2rad(lon[2] - lon[1])
   mflux(ii, jj, s) = rho[ii,jj,s,indt] * (zl[ii,jj,s,indt] - zl[ii,jj,s-1,indt])
   Threads.@threads :static for j in 1:Nj
      j0 = max(j - 1, 1); j1 = min(j + 1, Nj)
      cl = cosd(lat[j]); dy = Float32(R * deg2rad(lat[j1] - lat[j0]))
      Dc = Vector{Float32}(undef, Nk); mc = Vector{Float32}(undef, Nk)
      @inbounds for i in 1:Ni
         i0 = max(i - 1, 1); i1 = min(i + 1, Ni); dx = Float32(R * cl * dl * (i1 - i0))
         Dt = 0.0f0; Mt = 0.0f0
         for s in 2:Nk                                       # every layer, the top one included
            mx1 = mflux(i1, j, s) * U[i1,j,s,indt]; mx0 = mflux(i0, j, s) * U[i0,j,s,indt]
            my1 = mflux(i, j1, s) * V[i,j1,s,indt] * Float32(cosd(lat[j1]))
            my0 = mflux(i, j0, s) * V[i,j0,s,indt] * Float32(cosd(lat[j0]))
            D = (mx1 - mx0) / max(dx, 1.0f0) + (my1 - my0) / (max(dy, 1.0f0) * Float32(max(cl, 1.0e-3)))
            m = mflux(i, j, s)
            (isfinite(D) && isfinite(m)) || (D = 0.0f0; m = 0.0f0)
            Dc[s] = D; mc[s] = max(m, 0.0f0); Dt += D; Mt += mc[s]
         end
         lev[i,j,1,indt] = 0.0f0
         F = 0.0f0; Mb = 0.0f0
         for s in 2:Nk-1
            F -= Dc[s]; Mb += mc[s]
            Fs = F + (Mt > 0.0f0 ? Mb / Mt : 0.0f0) * Dt
            zfs  = 0.5f0 * (zl[i,j,s-1,indt] + zl[i,j,s,indt]); zfs1 = 0.5f0 * (zl[i,j,s,indt] + zl[i,j,s+1,indt])
            rh   = 0.5f0 * (rho[i,j,s,indt] + rho[i,j,s+1,indt])
            r = Fs / max(rh * (zfs1 - zfs), 1.0f-6)
            lev[i,j,s,indt] = isfinite(r) ? r : 0.0f0
         end
         lev[i,j,Nk,indt] = 0.0f0
      end
   end
   return nothing
end

# ── AROME ATMO-HDF (GRIB2, fixed height-above-ground levels, BL only) ────────
function metaromehdffieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   get(metconfig, "_layout_logged", false) || (metconfig["_layout_logged"] = true;
      @info "met field positions (AROME-HDF GRIB): every field on the regular grid at the same points (no staggering); " *
            "the fixed heights (and the stitched isobaric levels) are the full levels (slots 2…Nk, slot 1 = surface values); " *
            "zlev = their midpoints (half levels); vertical velocity: $(get(metconfig, "vertical_velocity", "continuity") == "continuity" ? "diagnosed from mass continuity on the half levels" : "none (0)")")
   Ni = Int(metconfig["Ni"]) + (globalx ? 1 : 0)
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   heights = Float32.(metconfig["levelvalues"])           # ascending [m AGL] (HP1 height block)
   stitch  = haskey(metconfig, "ip_pressures")            # isobaric package stitched on top?
   n_hp    = Int(get(metconfig, "n_hp", Nk - 1))          # # height levels (full column if no stitch)
   nhk     = n_hp + 1                                     # last height-block index (incl. surface slot)

   # Shift ONLY the height block up one to open the surface slot; the isobaric
   # block (slots nhk+1 … Nk, already read in place) must not be shifted.
   _shift_up_one!(metfields3D, ("Ux","Vy","T","rh","puv"), nhk, indt)

   zl = metfields3D["zlev"]
   zl[:,:,1,indt] .= 0.0f0
   for k in 2:nhk
      zl[:,:,k,indt] .= heights[k-1]
   end

   if haskey(metfields2D, "2t")
      metfields3D["T"][:,:,1,indt] = metfields2D["2t"][:,:,indt]
   else
      metfields3D["T"][:,:,1,indt] = metfields3D["T"][:,:,2,indt]
   end
   metfields3D["Ux"][:,:,1,indt] = metfields3D["Ux"][:,:,2,indt]
   metfields3D["Vy"][:,:,1,indt] = metfields3D["Vy"][:,:,2,indt]
   # Surface pressure ≈ lowest-level pressure barometrically extended to z=0.
   metfields3D["puv"][:,:,1,indt] = metfields3D["puv"][:,:,2,indt] .+
       _G0 .* 1.15f0 .* heights[1]

   # ── Pressure-level stitch (AROME-HDF IP1) ────────────────────────────────
   # Place the isobaric block at its geopotential height above ground.  Terrain
   # z_s from the hypsometric relation between MSL and surface pressure; the
   # ladder is forced strictly increasing across the seam.
   if stitch
      ippres = metconfig["ip_pressures"]::Vector{Float64}
      gz = metfields3D["gz"]                                 # geopotential [m² s⁻²] (IP1 z)
      T2 = haskey(metfields2D, "2t") ? metfields2D["2t"][:,:,indt] :
                                       metfields3D["T"][:,:,2,indt]
      prmsl = haskey(metfields2D, "prmsl") ? metfields2D["prmsl"][:,:,indt] : nothing
      sp    = haskey(metfields2D, "sp")    ? metfields2D["sp"][:,:,indt]    : nothing
      zs = (prmsl !== nothing && sp !== nothing) ?
           (287.05f0 .* max.(T2, 150.0f0) ./ _G0) .* log.(max.(prmsl, 1.0f0) ./ max.(sp, 1.0f0)) :
           zeros(Float32, Ni, Nj)
      for (r, p) in enumerate(ippres)
         k = nhk + r
         metfields3D["puv"][:,:,k,indt] .= Float32(p * 100.0)                  # hPa → Pa
         zl[:,:,k,indt] .= gz[:,:,k,indt] ./ _G0 .- zs                         # geopotential height AGL
      end
      # Isobaric levels at or below the top of the height block (700 hPa over flat
      # ground, down to ~450 hPa over the Alps) hold air the height block already
      # covers.  They used to be pushed 10 m above it with their own, lower-air
      # values: pressure rose with height in ~1.4 % of the ladder steps (2026-10-21,
      # check_met_readers).  Per column those slots now form a ramp between the top
      # height level and the first isobaric level above it.
      _arome_hdf_stitch_ramp!(metfields3D, nhk, Nk, Ni, Nj, indt)
   end

   # qv from relative humidity [%] over the FULL column (puv now set everywhere).
   if haskey(metfields3D, "rh") && haskey(metfields3D, "qv")
      T  = metfields3D["T"]
      p  = metfields3D["puv"]
      rh = metfields3D["rh"]
      tc = T[:,:,:,indt] .- 273.15f0
      es = 611.2f0 .* exp.(17.62f0 .* tc ./ (243.12f0 .+ tc))
      e  = clamp.(rh[:,:,:,indt], 0.0f0, 100.0f0) ./ 100.0f0 .* es
      metfields3D["qv"][:,:,:,indt] = 0.622f0 .* e ./ max.(p[:,:,:,indt] .- 0.378f0 .* e, 1.0f0)
   end

   # STACKED layout (2026-10-22): the fixed heights (and the stitched isobaric
   # levels) carry u, v, T, p — they are full levels; zlev goes to the half
   # levels between them (the transport reads zlev on half levels).
   _stack_half_levels!(zl, Ni, Nj, Nk, indt)
   _thermo_from_pTq!(metfields3D, Nk, indt)

   # Neither HP1 nor IP1 disseminates vertical velocity: diagnosed from mass
   # continuity (#info vertical_velocity = continuity, the default since
   # 2026-10-22) or zero (= none: vertical transport by turbulence only, the
   # behaviour before).  The isobaric stitch extends the horizontal winds and
   # thermodynamics to ~16 km.
   if get(metconfig, "vertical_velocity", "continuity") == "continuity"
      _continuity_levdot!(metconfig, metfields3D["levdot"], metfields3D["zlev"], metfields3D["rho"],
                          metfields3D["Ux"], metfields3D["Vy"], Ni, Nj, Nk, indt)
   else
      fill!(view(metfields3D["levdot"], :, :, :, indt), 0.0f0)
   end

   # SP3 sshf follows the ECMWF sign convention (positive downward); the
   # Obukhov-length derivation expects positive-upward flux.
   if haskey(metfields2D, "shf")
      metfields2D["shf"][:,:,indt] = -metfields2D["shf"][:,:,indt]
   end

   _surface_turbulence!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   return nothing
end

# ── MesoNH (NetCDF, Gal-Chen terrain-following levels, native TKE) ───────────
function metmesonhfieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   Ni = Int(metconfig["Ni"])
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   zhat = Float32.(metconfig["zhat"])                     # level heights [m], flat terrain
   ztop = zhat[end]

   # MesoNH is an Arakawa C grid (grid codes in the file: THT/PABST/TKET/RVT on
   # mass points = 1, UT on u-points = 2, VT on v-points = 3, ZHAT/WT on w-points
   # = 4).  ZHAT are the W-LEVEL heights (−10 m ghost, 0 m ground, …): mass
   # level k lies halfway between ZHAT[k] and ZHAT[k+1].  Until 2026-10-13 mass
   # level k was put at ZHAT[k] (the first physical level at the ground, clamped
   # to 1 m, instead of 5 m; hydrostatic balance off by +4.7 % median on the
   # 6 %-stretched AEROMARINE grid), WT was used as if on the mass levels, and
   # UT/VT on their faces as mass-point values (half a 2-km cell).
   # Gal-Chen height above ground: z(k,i,j) = ZHAT_mass[k]·(1 − ZS/ztop).
   # MesoNH level 1 is the below-ground ghost level → becomes the surface slot.
   zs = haskey(metfields2D, "zs") ? metfields2D["zs"][:,:,indt] : zeros(Float32, Ni, Nj)
   shrink = 1.0f0 .- zs ./ ztop
   # STACKED layout (2026-10-22; see _terrain_following_levdot_stacked!): zlev and
   # the vertical motion on the HALF levels = the w-levels, the mass-level fields
   # (u, v, θ, p, …, file level s = slot s) on the full levels between them.  The
   # full level of slot s (mass level s, between w-levels s and s+1) lies between
   # zlev[s-1] and zlev[s], so zlev[s] = ZHAT[s+1]: zlev[1] = ZHAT[2] = 0, the
   # ground (w-level 1 is the −10 m ghost).  2026-10-13 had put zlev on the mass
   # levels themselves, which the transport reads half a layer too low.
   zl = metfields3D["zlev"]
   zl[:,:,1,indt] .= 0.0f0
   for s in 2:Nk
      zh = s < Nk ? zhat[s+1] : zhat[Nk] + (zhat[Nk] - zhat[Nk-1])
      zl[:,:,s,indt] = max.(zh .* shrink, zl[:,:,s-1,indt] .+ 1.0f0)
   end
   # u, v from the faces to the mass points (u-point i = west face of cell i,
   # same count as the mass points; the last column/row keeps its face value)
   # face winds kept for the slope term of the coordinate velocity (C-grid form)
   # positions read from the file at setup (FieldPos: MesoNH `grid` 2/3 = faces,
   # 1 = mass points; 4 / level_w = w-levels): a post-processed file can carry
   # winds at the mass points or w on the mass levels
   fp = get(metconfig, "field_pos", Dict{String,FieldPos}())
   ux_face = !haskey(fp, "Ux") || fp["Ux"].x == :face
   vy_face = !haskey(fp, "Vy") || fp["Vy"].y == :face
   wz = haskey(fp, "levdot") ? fp["levdot"].z : :half
   tz = haskey(fp, "tke") ? fp["tke"].z : :full
   Uf = metfields3D["Ux"][:, :, :, indt:indt]; Vf = metfields3D["Vy"][:, :, :, indt:indt]
   for (key, dim) in (("Ux", 1), ("Vy", 2))
      haskey(metfields3D, key) || continue
      (dim == 1 ? ux_face : vy_face) || continue           # already at the mass points
      A = metfields3D[key]
      if dim == 1
         @views A[1:Ni-1,:,:,indt] .= 0.5f0 .* (A[1:Ni-1,:,:,indt] .+ A[2:Ni,:,:,indt])
      else
         @views A[:,1:Nj-1,:,indt] .= 0.5f0 .* (A[:,1:Nj-1,:,indt] .+ A[:,2:Nj,:,indt])
      end
   end

   # T from θ and absolute pressure; rho, drhodz.
   th = metfields3D["theta"]
   p  = metfields3D["puv"]
   metfields3D["T"][:,:,:,indt] = th[:,:,:,indt] .*
       (max.(p[:,:,:,indt], 1.0f0) ./ 100000.0f0).^0.28573f0
   # rv (mixing ratio) → qv (specific humidity)
   if haskey(metfields3D, "rv") && haskey(metfields3D, "qv")
      rv = max.(metfields3D["rv"][:,:,:,indt], 0.0f0)
      metfields3D["qv"][:,:,:,indt] = rv ./ (1.0f0 .+ rv)
   end
   tv = metfields3D["T"][:,:,:,indt]
   if haskey(metfields3D, "qv")
      tv = tv .* (1.0f0 .+ 0.608f0 .* metfields3D["qv"][:,:,:,indt])
   end
   metfields3D["rho"][:,:,:,indt] = p[:,:,:,indt] ./ (_RD .* max.(tv, 100.0f0))
   rho = metfields3D["rho"]
   for k in 2:Nk-1
      metfields3D["drhodz"][:,:,k,indt] = (rho[:,:,k+1,indt] .- rho[:,:,k-1,indt]) ./
          max.(zl[:,:,k+1,indt] .- zl[:,:,k-1,indt], 1.0f0)
   end
   metfields3D["drhodz"][:,:,1,indt]  = metfields3D["drhodz"][:,:,2,indt]
   metfields3D["drhodz"][:,:,Nk,indt] = metfields3D["drhodz"][:,:,Nk-1,indt]
   fill!(view(metfields3D["dsigwdz"], :, :, :, indt), 0.0f0)

   # WT [m/s] on the w-levels (w-level k = bottom of mass level k) → mass level
   # (ascending k: w[k+1] is still a w-level value when read) → levdot
   # [level-index/s] on the local dz ladder.
   w = metfields3D["levdot"]
   for s in 1:Nk-1                      # half level s = w-level s+1 (between mass levels s and s+1)
      if wz == :half
         @views w[:,:,s,indt] .= w[:,:,s+1,indt]
      else                              # w on the mass levels: mean of the two around it
         @views w[:,:,s,indt] .= 0.5f0 .* (w[:,:,s,indt] .+ w[:,:,s+1,indt])
      end
   end
   # coordinate velocity on the Gal-Chen half levels; 0 at the ground and the top
   # (C-grid slope term with the face winds; centred slopes when the file has
   # the winds at the mass points)
   if ux_face && vy_face
      wv = view(w, :, :, :, indt:indt); zv = view(zl, :, :, :, indt:indt)
      _terrain_following_levdot_stacked_cgrid!(metconfig, wv, zv, zs, Uf, Vf, Ni, Nj, Nk, 1)
   else
      _terrain_following_levdot_stacked!(metconfig, w, zl, zs, metfields3D["Ux"], metfields3D["Vy"], Ni, Nj, Nk, indt)
   end

   # MesoNH archives TKE directly (TKET → tke): the field is the FULL native
   # e = ½(σ_u²+σ_v²+σ_w²), NOT the Hanna 0.5·σ_w² fill — transport maps it via
   # the native-e closure (σ_w² = 2·e·r_w, K_z = C_K·ℓ·√e; MESONH.vars #info,
   # Fix C).  The Thomson scheme computes the turbulent length scale online, so
   # tls stays zero.
   if haskey(metfields3D, "tke")
      tk = metfields3D["tke"]
      replace!(view(tk, :, :, :, indt), NaN32 => 0.0f0)
      # TKET is on the mass levels (slot s = mass level s); the transport reads
      # tke on the half levels with zlev (stacked layout) → mean of the two mass
      # levels around each w-level (ascending: tk[s+1] is still a mass value)
      for s in 1:Nk-1
         if tz == :full
            @views tk[:,:,s,indt] .= 0.5f0 .* (tk[:,:,s,indt] .+ tk[:,:,s+1,indt])
         else                           # TKE on the w-levels: w-level s+1 = half level s
            @views tk[:,:,s,indt] .= tk[:,:,s+1,indt]
         end
      end
   end
   fill!(view(metfields3D["tls"], :, :, :, indt), 0.0f0)

   # Surface parameters (log-law ustar, Rib-based ol, TKE-threshold pblh) are
   # still derived — the Hanna option-b path and the PBL entrainment zone use
   # them — but the native TKE field is NOT overwritten by the Hanna fill.
   _derive_surface_params!(metfields3D, metfields2D, indt, Ni, Nj, Nk; native_tke=true)
   return nothing
end
