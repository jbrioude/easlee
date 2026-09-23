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
            z1  = max(zl[i,j,2,indt], 2.0f0)
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
      # Surface-layer bulk Richardson estimate (Monin-Obukhov similarity):
      # ζ ≈ Rib (unstable) or Rib/(1−5·Rib) (stable), evaluated at z(k=3).
      th = metfields3D["theta"]; u = metfields3D["Ux"]; v = metfields3D["Vy"]
      zl = metfields3D["zlev"]
      Threads.@threads :static for j in 1:Nj
         @inbounds for i in 1:Ni
            z   = zl[i,j,3,indt]
            dth = th[i,j,3,indt] - th[i,j,2,indt]
            du2 = (u[i,j,3,indt] - u[i,j,2,indt])^2 + (v[i,j,3,indt] - v[i,j,2,indt])^2 + 0.01f0
            (isfinite(dth) && isfinite(du2) && z > 1.0f0) || continue
            rib = _G0 * z * dth / (max(th[i,j,2,indt], 200.0f0) * du2)
            zeta = rib >= 0.0f0 ? rib / max(1.0f0 - 5.0f0 * rib, 0.1f0) : rib
            metfields2D["ol"][i,j,indt] = abs(zeta) < 1.0f-4 ? 5000.0f0 :
                                          clamp(z / zeta, -5000.0f0, 5000.0f0)
         end
      end
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
# lowest height where Rib(z) = g·z·(θv(z)−θv_sfc) / (θv_sfc·|ΔU|²) crosses Ri_c.
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
         ribp = 0.0f0; zp = zl[i,j,2,indt]
         for k in 3:Nk
            z = zl[i,j,k,indt]
            thv = th[i,j,k,indt] * (qv === nothing ? 1.0f0 : 1.0f0 + 0.61f0 * max(qv[i,j,k,indt], 0.0f0))
            du2 = (u[i,j,k,indt] - u[i,j,2,indt])^2 + (v[i,j,k,indt] - v[i,j,2,indt])^2 + 0.1f0
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

   _thermo_from_pTq!(metfields3D, Nk, indt)

   # OMEGA [Pa/s] → levdot [level-index/s]:  dindl/dp = 2/(p[k+1]−p[k−1]) < 0
   p = metfields3D["puv"]
   for k in 2:Nk-1
      dp = p[:,:,k+1,indt] .- p[:,:,k-1,indt]
      metfields3D["levdot"][:,:,k,indt] = metfields3D["levdot"][:,:,k,indt] .* 2.0f0 ./
          min.(dp, -1.0f0)
   end
   metfields3D["levdot"][:,:,1,indt]  = metfields3D["levdot"][:,:,2,indt] .* 0.5f0
   metfields3D["levdot"][:,:,Nk,indt] .= 0.0f0

   # Total precipitation flux for wet deposition: sum of the mst_Ne edge fluxes
   # (convective + large-scale/anvil, rain + ice).  Edge vars are NOT shifted.
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

# ── NCEP GFS (GRIB2, isobaric levels) ─────────────────────────────────────────
function metgfsfieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
   Ni = Int(metconfig["Ni"]) + (globalx ? 1 : 0)
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   plev = Float32.(metconfig["levelvalues"]) .* 100.0f0   # hPa → Pa, descending

   _shift_up_one!(metfields3D, ("Ux","Vy","T","qv","levdot","zlev","cc"), Nk, indt)

   sp   = metfields2D["sp"][:,:,indt]
   orog = haskey(metfields2D, "orog") ? metfields2D["orog"][:,:,indt] : zeros(Float32, Ni, Nj)

   # Heights above ground from geopotential height; enforce monotonicity where
   # pressure levels run below the terrain.
   zl = metfields3D["zlev"]
   zl[:,:,1,indt] .= 0.0f0
   for k in 2:Nk
      zl[:,:,k,indt] = zl[:,:,k,indt] .- orog
      zl[:,:,k,indt] = max.(zl[:,:,k,indt], zl[:,:,k-1,indt] .+ 1.0f0)
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

   _thermo_from_pTq!(metfields3D, Nk, indt)

   # omega [Pa/s] → levdot, using the FIXED isobaric ladder spacing (never the
   # surface pressure at k=1: over high terrain sp can drop below the first
   # levels and the centred dp would degenerate).
   for k in 2:Nk-1
      dp = k == 2 ? 2.0f0 * (plev[2] - plev[1]) : plev[k] - plev[k-2]   # < 0
      metfields3D["levdot"][:,:,k,indt] = metfields3D["levdot"][:,:,k,indt] .* (2.0f0 / dp)
   end
   metfields3D["levdot"][:,:,1,indt]  = metfields3D["levdot"][:,:,2,indt] .* 0.5f0
   metfields3D["levdot"][:,:,Nk,indt] .= 0.0f0

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

# ── AROME ATMO-HDF (GRIB2, fixed height-above-ground levels, BL only) ────────
function metaromehdffieldcalculation!(globalx::Bool, metconfig, metfields3D, metfields2D, indt)
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
         hk = gz[:,:,k,indt] ./ _G0 .- zs                                      # geopotential height AGL
         zl[:,:,k,indt] = max.(hk, zl[:,:,k-1,indt] .+ 10.0f0)                 # keep monotonic
      end
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

   _thermo_from_pTq!(metfields3D, Nk, indt)

   # Neither HP1 nor IP1 disseminates vertical velocity → levdot stays zero;
   # vertical transport is turbulence-only.  The isobaric stitch extends the
   # horizontal winds + thermodynamics to ~16 km (no longer a 3 km lid).
   fill!(view(metfields3D["levdot"], :, :, :, indt), 0.0f0)

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

   # Gal-Chen height above ground: z(k,i,j) = ZHAT[k]·(1 − ZS/ztop).
   # MesoNH level 1 is the below-ground ghost level → becomes the surface slot.
   zs = haskey(metfields2D, "zs") ? metfields2D["zs"][:,:,indt] : zeros(Float32, Ni, Nj)
   shrink = 1.0f0 .- zs ./ ztop
   zl = metfields3D["zlev"]
   zl[:,:,1,indt] .= 0.0f0
   for k in 2:Nk
      zl[:,:,k,indt] = max.(zhat[k] .* shrink, zl[:,:,k-1,indt] .+ 1.0f0)
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

   # WT [m/s] → levdot [level-index/s] on the local dz ladder.
   w = metfields3D["levdot"]
   for k in 2:Nk-1
      dz = (zl[:,:,k+1,indt] .- zl[:,:,k-1,indt]) ./ 2.0f0
      w[:,:,k,indt] = w[:,:,k,indt] ./ max.(dz, 1.0f0)
   end
   w[:,:,1,indt]  = w[:,:,2,indt] .* 0.5f0
   w[:,:,Nk,indt] .= 0.0f0

   # MesoNH archives TKE directly (TKET → tke): the field is the FULL native
   # e = ½(σ_u²+σ_v²+σ_w²), NOT the Hanna 0.5·σ_w² fill — transport maps it via
   # the native-e closure (σ_w² = 2·e·r_w, K_z = C_K·ℓ·√e; MESONH.vars #info,
   # Fix C).  The Thomson scheme computes the turbulent length scale online, so
   # tls stays zero.
   if haskey(metfields3D, "tke")
      replace!(view(metfields3D["tke"], :, :, :, indt), NaN32 => 0.0f0)
   end
   fill!(view(metfields3D["tls"], :, :, :, indt), 0.0f0)

   # Surface parameters (log-law ustar, Rib-based ol, TKE-threshold pblh) are
   # still derived — the Hanna option-b path and the PBL entrainment zone use
   # them — but the native TKE field is NOT overwritten by the Hanna fill.
   _derive_surface_params!(metfields3D, metfields2D, indt, Ni, Nj, Nk; native_tke=true)
   return nothing
end
