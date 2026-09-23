# ──────────────────────────────────────────────────────────────────────────────
# io/initialize.jl — met-input ingestion & runtime initialisation
#
# Scans and loads the meteorological input (GRIB via eccodes; NetCDF for
# WRF / AROME / MERRA2 / MesoNH), builds the metconfig / gridconfig /
# levelconfig dictionaries, allocates the two-time-slot met arrays, and
# derives the fields the transport kernels need (hydrostatic zlev/puv
# integration, rho/theta/levdot, Hanna surface-turbulence parameters).
# FLEXPART lineage: readwind / readavailable / verttransform / gridcheck.
#
# Entry points (called from runtime.jl, build_run.jl and prefetch.jl):
#   gridcheck_onefile!   — inspect one reference met file → metconfig dict
#   initializerun / initoutgrid / initgridconfig / initlevelconfig / initarrays
#   readmetfiles!        — keep the two met time slots in step with the run
#   readflexpartinput    — legacy FLEXPART option-file ingestion (build_run.jl)
#
# Includes at the bottom: metcalc_models.jl, vmix_namelist.jl, prefetch.jl.
# ──────────────────────────────────────────────────────────────────────────────

using GRIB
using Dates
using Interpolations
using Statistics

function netcdfinfo_WRF(filename)
   info = Dict()
   print("filename netcdf=",filename)
   v = NetCDF.open(filename)
   info["dy"] = v.gatts["DY"] #in meters
   info["dx"] = v.gatts["DX"] #in meters
   info["Ni"] = v.gatts["WEST-EAST_PATCH_END_UNSTAG"]
   info["Nj"] = v.gatts["SOUTH-NORTH_PATCH_END_UNSTAG"]
   info["Nk"] = v.gatts["BOTTOM-TOP_PATCH_END_STAG"] #for U,V,T, the first level will be the surface anyway
   edge = v["U"].atts["stagger"]
   if edge == "X"
      info["Uedge"] = true
   else
      info["Uedge"] = false
   end
   edge = v["V"].atts["stagger"]
   if edge == "Y"
      info["Vedge"] = true
   else
      info["Vedge"] = false
   end
   edge = v["W"].atts["stagger"]
   if edge == "Z"
      info["Zedge"] = true
   else
      info["Zedge"] = false
   end
   info["CEN_LAT"] = v.gatts["CEN_LAT"]
   info["CEN_LON"] = v.gatts["CEN_LON"]
   info["proj_truelat1"] = v.gatts["TRUELAT1"]
   info["proj_truelat2"] = v.gatts["TRUELAT2"]
   info["GRID_ID"] = v.gatts["GRID_ID"]
   info["proj_stdlon"] = v.gatts["STAND_LON"]
   info["I_PARENT_START"] = v.gatts["I_PARENT_START"]
   info["J_PARENT_START"] = v.gatts["J_PARENT_START"]
   info["PARENT_ID"] = v.gatts["PARENT_ID"]
   info["MAP_PROJ"] = v.gatts["MAP_PROJ"]
   if info["MAP_PROJ"] == 1
      info["projection"] = "Lambert"
   elseif info["MAP_PROJ"] == 2
      info["projection"] = "Polarstereo"
   elseif info["MAP_PROJ"] == 3
      info["projection"] = "mercator"
   elseif info["MAP_PROJ"] == 4
      info["projection"] = "Rectilinear"
   end 
   proj_truelat1 = info["proj_truelat1"]
   info["longitudes"] = v.vars["XLONG"][:,:,1]
   info["latitudes"] = v.vars["XLAT"][:,:,1]
   proj_lat1 = info["latitudes"][1,1]
   proj_lon1 = info["longitudes"][1,1]
   proj_dx = info["dx"]
   earth_radius_m = 6370000
   proj_rebydx = earth_radius_m / proj_dx
   if (info["proj_truelat1"] < 0.)
      proj_hemi = -1.0
   else
      proj_hemi = 1.0
   end
   proj_cone, proj_rsw, proj_polei, proj_polej = set_lc(info["proj_truelat1"],info["proj_truelat2"], proj_lon1, info["proj_stdlon"], proj_rebydx, proj_hemi, proj_lat1)
   info["proj_hemi"] = proj_hemi
   info["proj_cone"] = proj_cone
   info["proj_rsw"] = proj_rsw
   info["proj_polei"] = proj_polei
   info["proj_polej"] = proj_polej
   info["proj_rebydx"] = proj_rebydx
   info["proj_lat1"] = proj_lat1
   info["proj_rebydx"] = proj_rebydx
   return info 
end
function netcdfinfo_AROME(filename)
   info = Dict()
   print("filename netcdf=",filename)
   v = NetCDF.open(filename)
   info["dy"] = v.gatts["y_resolution"] #in meters
   info["dx"] = v.gatts["x_resolution"] #in meters
#akbk from AROME
   info["ak"] = v.vars["Ak"]
   info["bk"] = v.vars["Bk"]
   # Grid sizes are hardcoded for this AROME domain (the files lack the
   # WRF-style dimension attributes).  Nk = file levels + 1 because of NH_dW.
   info["Ni"] = 220
   info["Nj"] = 190
   info["Nk"] = 91
   # AROME files carry no stagger attributes (cf. the WRF branch): U/V live on
   # mass points, W on level edges.
   info["Uedge"] = false
   info["Vedge"] = false
   info["Zedge"] = true
   info["proj_truelat1"] = v.gatts["center_lat"]
   info["proj_truelat2"] = v.gatts["center_lat"]
   info["proj_stdlon"] = v.gatts["reference_lon"]
   info["MAP_PROJ"] = 1   # AROME native grid: Lambert conformal
   if info["MAP_PROJ"] == 1
      info["projection"] = "Lambert"
   elseif info["MAP_PROJ"] == 2
      info["projection"] = "Polarstereo"
   elseif info["MAP_PROJ"] == 3
      info["projection"] = "mercator"
   elseif info["MAP_PROJ"] == 4
      info["projection"] = "Rectilinear"
   end 
   proj_truelat1 = info["proj_truelat1"]
   info["longitudes"] = v.vars["LON"][:,:,1]
   info["latitudes"] = v.vars["LAT"][:,:,1]
   proj_lat1 = info["latitudes"][1,1]
   proj_lon1 = info["longitudes"][1,1]
   proj_dx = info["dx"]
   earth_radius_m = 6371229
   proj_rebydx = earth_radius_m / proj_dx
   if (info["proj_truelat1"] < 0.)
      proj_hemi = -1.0
   else
      proj_hemi = 1.0
   end
   proj_cone, proj_rsw, proj_polei, proj_polej = set_lc(info["proj_truelat1"],info["proj_truelat2"], proj_lon1, info["proj_stdlon"], proj_rebydx, proj_hemi, proj_lat1)
   info["proj_hemi"] = proj_hemi
   info["proj_cone"] = proj_cone
   info["proj_rsw"] = proj_rsw
   info["proj_polei"] = proj_polei
   info["proj_polej"] = proj_polej
   info["proj_rebydx"] = proj_rebydx
   info["proj_lat1"] = proj_lat1
   info["proj_rebydx"] = proj_rebydx
   return info 
end
function netcdfinfo_MERRA2(filename)
   info = Dict()
   lons = Float32.(NetCDF.ncread(filename, "lon"))
   lats = Float32.(NetCDF.ncread(filename, "lat"))
   lev  = NetCDF.ncread(filename, "lev")
   info["Ni"] = Int32(length(lons))
   info["Nj"] = Int32(length(lats))
   # +1 surface slot; the mst_Ne edge collection (73 entries) fits exactly.
   info["Nk"] = Int32(length(lev) + 1)
   info["longitudes"] = lons
   info["latitudes"]  = lats
   info["dx"] = lons[2] - lons[1]
   info["dy"] = lats[2] - lats[1]
   info["projection"] = "Rectilinear"
   dlon = lons[2] - lons[1]
   span = mod(maximum(lons) - minimum(lons), 360.0f0)
   info["globalx"] = span == 0.0f0 || isapprox(span, 360.0f0 - dlon; atol=1e-3)
   return info
end

function netcdfinfo_MESONH(filename)
   info = Dict()
   xhat = NetCDF.ncread(filename, "XHAT")    # [m] u-point x coordinates
   yhat = NetCDF.ncread(filename, "YHAT")
   zhat = NetCDF.ncread(filename, "ZHAT")    # [m] w-level heights, flat terrain
   lat2d = Float32.(NetCDF.ncread(filename, "LAT"))
   lon2d = Float32.(NetCDF.ncread(filename, "LON"))
   rpk  = NetCDF.ncread(filename, "RPK")[1]
   beta = NetCDF.ncread(filename, "BETA")[1]
   abs(beta) < 1e-6 || error("MESONH: rotated conformal projection (BETA=$beta) not supported")
   abs(rpk)  < 1e-6 || error("MESONH: Lambert conformal (RPK=$rpk) not wired yet — only Mercator (RPK=0)")
   info["Ni"] = Int32(size(lat2d, 1))
   info["Nj"] = Int32(size(lat2d, 2))
   info["Nk"] = Int32(length(zhat))          # MesoNH ghost level 1 = surface slot
   info["zhat"] = zhat
   info["projection"] = "Mercator"
   info["dx"] = Float32(xhat[2] - xhat[1])
   info["dy"] = Float32(yhat[2] - yhat[1])
   info["proj_lat0"]  = Float32(NetCDF.ncread(filename, "LAT0")[1])  # reference latitude
   info["proj_lator"] = Float32(lat2d[1, 1]) # first mass point
   info["proj_lonor"] = Float32(lon2d[1, 1])
   info["longitudes"] = lon2d
   info["latitudes"]  = lat2d
   info["globalx"] = false
   return info
end

function netcdfgetdims(filename,input_vars_3D,input_vars_2D)
   dims = Dict()
   v = NetCDF.open(filename)
   for name in input_vars_3D
     dims[name] = size(v[name])[1:3]
   end
   for name in input_vars_2D
     dims[name] = size(v[name])[1:2]
   end
   return dims
end
"""
    gridcheck_onefile!(directory, input) -> metconfig::Dict

Inspect ONE reference met file (GRIB or NetCDF) and build the `metconfig`
dictionary: grid geometry and projection, vertical-coordinate ladder (ak/bk or
level values), variable-name mappings from the variables namelist, companion
file patterns, and the vertical-mixing closure.  FLEXPART lineage: gridcheck.
"""
function gridcheck_onefile!(directory,input)
   list_metconfig = readdir(directory,join=true)
   list_metconfig = filter(x -> !occursin(r"AVAILABLE", x), list_metconfig)
   local metconfig = Dict()
   metconfig["filename"] = input["filename"] #list of filenames now in metconfig
   metconfig["filedate"] = input["filedate"] #list of file dates now in metconfig
   metconfig["leveltype"] = nothing 
   metconfig["center"] = nothing 
   metconfig["gridtype"] = nothing 
   metconfig["metfileformat"] = nothing 
   metconfig["globalx"] = false
   convertnameinput = nothing
   convertnameeaslee = nothing
   easlee_vars_3D = nothing
   input_vars_3D = nothing
   easlee_vars_2D = nothing
   input_vars_2D = nothing
   # Grid reference file: prefer the first AVAILABLE entry (robust for nested
   # directory layouts like GFS date/cycle subdirs, and guarantees we inspect a
   # main file rather than a companion).  Fall back to the directory listing.
   filename = isempty(input["filename"]) ? list_metconfig[1] : input["filename"][1]
   convertnameeaslee, convertnameinput, easlee_vars_3D, input_vars_3D, easlee_vars_2D, input_vars_2D = extract_variables!(input)
   metconfig["metfileformat"] = uppercase(strip(input["format"]))

   # ── Fix C: vertical-mixing closure (σ_w² / K_z semantics of the tke field) ─
   # Selected at read time from the namelist #info keys (tke_kind / vmix_rw /
   # vmix_ck — see vmix_namelist.jl); errors here, before any transport, when a
   # native-TKE driver is missing its closure constants (incl. the WRF EXCH_H
   # trap).  Published to the kernels by build_runtime.
   tke_mapped = strip(string(get(convertnameeaslee, "tke", ""))) != ""
   kz_mapped  = strip(string(get(convertnameeaslee, "kz", ""))) != "" ||
                strip(string(get(convertnameeaslee, "exch_h", ""))) != ""
   metconfig["vmix"] = vmix_from_namelist(input, tke_mapped, kz_mapped)

   # ── Vertical-velocity convention of the archived "levdot" field ───────────
   # TWO conventions exist for the SAME parameter (ECMWF 77):
   #   raw MARS / ERA5 dissemination : eta-dot            [s^-1]
   #   flex_extract (calc_etadot)    : eta-dot * dp/deta  [Pa s^-1]
   # They differ by dp/deta ~ 1e5, so the conversion to a level-index rate must
   # divide by the eta thickness in the first case and by dp in the second.
   # Declared per dataset with the #info key `etadot_units` (auto | s-1 | pa/s);
   # "auto" decides per file from the field magnitude, which is unambiguous at
   # five orders of separation.  See docs/ETADOT_UNITS.md.
   metconfig["etadot_units"] =
       lowercase(strip(string(get(input, "etadot_units", "auto"))))
   metconfig["etadot_units"] in ("auto", "s-1", "pa/s") ||
       error("listvariables #info etadot_units must be auto | s-1 | pa/s, got " *
             "\"$(metconfig["etadot_units"])\"")
   # Decided ONCE per met source: a file set is either flex_extract output or
   # raw dissemination, never a mix.  `nothing` = not yet determined (auto).
   metconfig["etadot_per_second"] =
       metconfig["etadot_units"] == "s-1"  ? true  :
       metconfig["etadot_units"] == "pa/s" ? false : nothing

   # ── Unit convention of the archived surface fluxes (sshf, ewss/nsss) ──────
   # review 2026-09-04 A3 (item 1): the GRIB header (stepType=accum, units
   # J m**-2 / N m**-2 s) cannot be trusted — flex_extract de-accumulates the
   # VALUES to per-second and keeps the accumulated header, raw MARS keeps both
   # accumulated.  Declared per dataset with `#info sfcflux_units`
   # (auto | per-s | accum); "auto" decides ONCE per met source from the field
   # magnitude (see _sfcflux_convention!), like etadot_units above.
   metconfig["sfcflux_units"] =
       lowercase(strip(string(get(input, "sfcflux_units", "auto"))))
   metconfig["sfcflux_units"] in ("auto", "per-s", "accum") ||
       error("listvariables #info sfcflux_units must be auto | per-s | accum, got " *
             "\"$(metconfig["sfcflux_units"])\"")
   metconfig["sfcflux_accum"] =
       metconfig["sfcflux_units"] == "accum" ? true  :
       metconfig["sfcflux_units"] == "per-s" ? false : nothing

   # ── Optional namelist #info keys (multi-source support) ───────────────────
   # companions: filename-pattern substitutions for datasets that split fields
   #             across several files at one valid time,
   #             e.g.  companions = _MET_>_CLD_
   #             or    companions = _asm_Nv_>_mst_Ne_, _asm_Nv_>_trb_Ne_
   # leveltype : vertical coordinate of the 3D GRIB messages
   #             ("hybrid" default, "isobaric", "height")
   # multistep : "true" when one GRIB file holds several valid times
   #             (messages are then filtered by validityDate/validityTime)
   comp = Tuple{String,String}[]
   if haskey(input, "companions")
      for pat in split(input["companions"], ",")
         parts = split(strip(pat), ">")
         length(parts) == 2 && push!(comp, (String(strip(parts[1])), String(strip(parts[2]))))
      end
   end
   metconfig["companions"] = comp
   # pressure_companion: a 3D isobaric package stitched ABOVE a height ladder
   # (AROME-HDF: HP1 height levels 0–3 km + IP1 pressure levels higher up).
   pc = nothing
   if haskey(input, "pressure_companion")
      parts = split(strip(input["pressure_companion"]), ">")
      length(parts) == 2 && (pc = (String(strip(parts[1])), String(strip(parts[2]))))
   end
   metconfig["pressure_companion"] = pc
   metconfig["vcoord"]     = lowercase(strip(get(input, "leveltype", "hybrid")))
   metconfig["multistep"]  = lowercase(strip(get(input, "multistep", "false"))) == "true"
   if metconfig["metfileformat"] == "GRIB"
     message = GribFile(filename)
     info = read(message)
     metconfig["projection"] = "Rectilinear"
     metconfig["gridtype"] = info["gridType"]
     metconfig["leveltype"] = info["typeOfLevel"]
     metconfig["center"] = info["typeOfLevel"]
     GRIB.destroy(message)
     metconfig = gribinfo(info,metconfig)
     if metconfig["vcoord"] != "hybrid"
        # Non-hybrid GRIB (isobaric or height levels): build the level-value
        # ladder by scanning the reference file.  Use the zonal-wind variable
        # only — other vars may carry extra levels (e.g. AROME-HDF t at 2 m)
        # that the dynamics fields don't have.
        tol = metconfig["vcoord"] == "isobaric" ? "isobaricInhPa" : "heightAboveGround"
        uxname = convertnameeaslee["Ux"]
        levvals = Float64[]
        GribFile(filename) do f
           for msg in f
              m = Message(msg)
              (m["shortName"] == uxname && m["typeOfLevel"] == tol) || continue
              lv = Float64(m["level"])
              lv in levvals || push!(levvals, lv)
           end
        end
        sort!(levvals, rev = (metconfig["vcoord"] == "isobaric"))
        metconfig["levelvalues"] = levvals
        metconfig["leveltol"]    = tol
        metconfig["Nk"] = Int32(length(levvals) + 1)   # +1 surface slot

        # Pressure-level stitch: append the isobaric package's levels (above the
        # ~3 km height cap → p ≲ 700 hPa) as extra ladder slots.  Their per-column
        # heights are filled later from the geopotential; here we only size Nk.
        if pc !== nothing && occursin(pc[1], basename(filename))
           ipfile = joinpath(dirname(filename), replace(basename(filename), pc[1] => pc[2]))
           if isfile(ipfile)
              ippres = Float64[]
              GribFile(ipfile) do f
                 for msg in f
                    m = Message(msg)
                    (m["shortName"] == uxname && m["typeOfLevel"] == "isobaricInhPa") || continue
                    p = Float64(m["level"])
                    (p <= 700.0 && !(p in ippres)) && push!(ippres, p)
                 end
              end
              sort!(ippres, rev = true)               # descending p = ascending height
              metconfig["ip_pressures"] = ippres
              metconfig["n_hp"]         = length(levvals)
              metconfig["Nk"]           = Int32(length(levvals) + 1 + length(ippres))
              @info "AROME-HDF pressure stitch: $(length(levvals)) height levels + $(length(ippres)) isobaric levels (≤700 hPa) → Nk=$(metconfig["Nk"])"
           end
        end
     end
   end
   if metconfig["metfileformat"] == "NETCDF"
     modeltype = uppercase(strip(input["model"]))
     # For WRF, dimensions are extracted from the NetCDF file using variable names specific to WRF output.
     # The info dictionary is merged to include grid and projection details, and inputdimensions are set for later array allocation.
     if modeltype == "WRF"
       info = netcdfinfo_WRF(filename)
       metconfig = merge(metconfig,info)
       dimensions = netcdfgetdims(filename,input_vars_3D,input_vars_2D)
       metconfig["inputdimensions"] = dimensions
     end  
     # For AROME, dimensions and grid details are handled differently due to the model's unique structure and variable naming.
     # The info dictionary includes ak/bk coefficients and fixed grid sizes, and inputdimensions are set for array allocation.
     if modeltype == "AROME"
       info = netcdfinfo_AROME(filename)
       metconfig = merge(metconfig,info)
       dimensions = netcdfgetdims(filename,input_vars_3D,input_vars_2D)
       metconfig["inputdimensions"] = dimensions
     end
     # MERRA2 and MESONH use the generic NetCDF path: uniform grids, no
     # inputdimensions entry (arrays allocated (Ni[+1], Nj, Nk, 2) directly).
     if modeltype == "MERRA2"
       info = netcdfinfo_MERRA2(filename)
       metconfig = merge(metconfig,info)
     end
     if modeltype == "MESONH"
       info = netcdfinfo_MESONH(filename)
       metconfig = merge(metconfig,info)
     end
   end
   metconfig["convertnameeaslee"] = convertnameeaslee
   metconfig["convertnameinput"] = convertnameinput
   metconfig["easlee_vars_3D"] = easlee_vars_3D
   metconfig["input_vars_3D"] = input_vars_3D
   metconfig["easlee_vars_2D"] = easlee_vars_2D
   metconfig["input_vars_2D"] = input_vars_2D
   metconfig["model"] = input["model"]
   return metconfig
end

function gribinfo(test,metconfig)
   if test["typeOfLevel"] == "hybrid"
      akbk = test["pv"]
      metconfig["Nk"] = Int32(length(akbk)/2)
      Nk = metconfig["Nk"]
      metconfig["ak"] = akbk[1:Nk]
      metconfig["bk"] = akbk[Nk+1:end]
   end
   metconfig["Ni"] = test["Ni"]
   metconfig["Nj"] = test["Nj"]
   metconfig["longitudes"] = reshape(test["longitudes"],test["Ni"],test["Nj"])
   # review 2026-09-04 A3 (item 5): flip latitude to ascending only when the scan
   # is north→south (jScansPositively == 0, the ECMWF/GFS default); a south→north
   # archive is stored as is.  The same flag drives the message store.
   flip_lat = true
   try; flip_lat = Int(test["jScansPositively"]) == 0; catch; end
   metconfig["flip_lat"] = flip_lat
   lat2 = reshape(test["latitudes"],test["Ni"],test["Nj"])
   metconfig["latitudes"] = flip_lat ? reverse(lat2,dims=2) : lat2
   metconfig["longitudes"] = metconfig["longitudes"][:,1]
   metconfig["latitudes"] = metconfig["latitudes"][1,:]
   issorted(metconfig["latitudes"]; lt = <) && allunique(metconfig["latitudes"]) ||
      error("GRIB latitude axis is not strictly ascending after the scan-direction flip " *
            "(jScansPositively=$(flip_lat ? 0 : 1), lat[1]=$(metconfig["latitudes"][1]), " *
            "lat[end]=$(metconfig["latitudes"][end]))")
   lon1 = minimum(metconfig["longitudes"])
   lon2 = maximum(metconfig["longitudes"]) 
   dlon = metconfig["longitudes"][2]-metconfig["longitudes"][1]
   if mod(lon2-lon1,360)==0 || mod(lon2-lon1,360) == dlon || mod(lon2-lon1,360) == 360-dlon
      metconfig["globalx"] = true
   else
      metconfig["globalx"] = false
   end

   return metconfig
end #end gribinfo

# ── Domain crop (domain-fill runs) ───────────────────────────────────────────
# Shrink the met grid to the release box + `halo` cells BEFORE the field arrays
# are allocated: overwrite the axes/Ni/Nj in metconfig and record the gather
# window; the GRIB reader then stores only the sub-window of each decoded field.
# Everything downstream (initgridconfig → Rectilinear origin/size, interp,
# advance, output mapping) is origin-relative, so the crop is transparent.
#
# The window may cross a global grid's longitude seam (Pacific box 150–230°E on
# a −180..180 ERA5 frame): the gather index list `iidx` rolls it contiguous and
# the new longitude axis continues monotonically past the seam (e.g. 150…229.75)
# — latlontoindyindx's mod-normalisation maps any parcel/receptor lon into that
# frame.  A lon-cropped grid is no longer zonally periodic → globalx = false
# (no wrap column, no periodic advance wrap; the domain-fill cull keeps parcels
# inside the box anyway).  The halo only needs to cover the 2-point
# interpolation stencil (1 cell); 2 leaves margin.  GRIB met only — for other
# formats the crop is skipped with a warning.
function _crop_metconfig!(metconfig, lon0, lat0, lon1, lat1; halo::Int = 2)
   if metconfig["metfileformat"] != "GRIB"
      @warn "[meteorology].crop: only implemented for GRIB met — running uncropped" format=metconfig["metfileformat"]
      return metconfig
   end
   lons = Float64.(metconfig["longitudes"]);  lats = Float64.(metconfig["latitudes"])
   Ni = Int(metconfig["Ni"]);  Nj = Int(metconfig["Nj"])
   dlon = lons[2] - lons[1];   dlat = lats[2] - lats[1]      # both ascending axes
   (dlon > 0 && dlat > 0) || error("crop: met axes not ascending (dlon=$dlon dlat=$dlat)")

   # latitude window (simple clamp — no wrap in latitude)
   la0 = min(lat0, lat1); la1 = max(lat0, lat1)
   (la0 <= lats[end] && la1 >= lats[1]) ||
      error("crop: latitude box [$la0, $la1] lies outside the met axis $(lats[1])…$(lats[end])")
   j0 = clamp(searchsortedlast(lats, la0) - halo, 1, Nj)
   j1 = clamp(searchsortedfirst(lats, la1) + halo, 1, Nj)
   j1 > j0 || error("crop: empty latitude window [$la0, $la1] on axis $(lats[1])…$(lats[end])")

   iidx::Union{Nothing,Vector{Int32}} = nothing
   local i0::Int
   if Bool(metconfig["globalx"])
      # periodic axis: unrolled indices may pass the seam (roll-to-contiguous)
      u0 = mod(lon0 - lons[1], 360.0) / dlon + 1.0         # west edge, fractional
      u1 = u0 + mod(lon1 - lon0, 360.0) / dlon             # east edge ≥ west edge
      i0 = floor(Int, u0) - halo
      i1 = ceil(Int, u1) + halo
      (i1 - i0 + 1) < Ni &&
         (iidx = Int32[mod1(i, Ni) for i in i0:i1])        # else: full circle — keep lon
   else
      # regional axis: normalise box lons into the grid-centred ±180° frame
      # (a −84° box edge on a 260–300°E grid → 276°E; no wrap past the edges)
      l0n = lons[1] + mod(lon0 - lons[1] + 180.0, 360.0) - 180.0
      l1n = lons[1] + mod(lon1 - lons[1] + 180.0, 360.0) - 180.0
      l1n < l0n && (l1n = l0n + mod(lon1 - lon0, 360.0))
      (l0n <= lons[end] && l1n >= lons[1]) ||
         error("crop: longitude box [$lon0, $lon1] lies outside the met axis $(lons[1])…$(lons[end])")
      i0 = clamp(floor(Int, (l0n - lons[1]) / dlon + 1.0) - halo, 1, Ni)
      i1 = clamp(ceil(Int,  (l1n - lons[1]) / dlon + 1.0) + halo, 1, Ni)
      (i1 - i0 + 1) < Ni && (iidx = collect(Int32, i0:i1))
   end

   if iidx !== nothing
      metconfig["longitudes"] = [lons[1] + (i - 1) * dlon for i in i0:(i0 + length(iidx) - 1)]
      metconfig["Ni"]      = Int32(length(iidx))
      metconfig["globalx"] = false                          # a sub-window is not periodic
   else
      iidx = collect(Int32, 1:Ni)                           # lat-only crop: gather all lons
   end
   metconfig["latitudes"] = lats[j0:j1]
   metconfig["Nj"]        = Int32(j1 - j0 + 1)
   metconfig["crop"]      = (iidx, j0:j1)
   @info "met crop: box ($lon0…$(lon1)E, $lat0…$(lat1)N) + $halo cells → $(metconfig["Ni"])×$(metconfig["Nj"]) of $(Ni)×$(Nj) ($(round(100 * Int(metconfig["Ni"]) * Int(metconfig["Nj"]) / (Ni * Nj), digits=1))% of the full grid)"
   return metconfig
end

function readnamelist(file::AbstractString)
    namelist_dict = Dict{String, Any}()
    keyold = ""
    open(file) do f
        current_namelist = ""
        for line in eachline(f)
            line = strip(line)
            if occursin(r"^&", line)
                current_namelist = match(r"^&(.*)", line).captures[1]
                if !haskey(namelist_dict, current_namelist)
                   namelist_dict[current_namelist] = Dict{String, Any}()
                end
            elseif occursin(r"^[a-zA-Z]", line)
                key, value = match(r"(.*)=(.*)", line).captures
                key = strip(key)
                value = strip(value[1:end-1])
                if !haskey(namelist_dict[current_namelist],key)
                    namelist_dict[current_namelist][key] = [value,]
                else
                  push!(namelist_dict[current_namelist][key],value)
                end
                keyold = key
            elseif occursin(r"^[0-9]",line)
                # Continuation line (bare numbers): append to the previous key.
                line = ","*line
                push!(namelist_dict[current_namelist][keyold],line)
            end
        end
    end
    return namelist_dict
end

function readflexpartinput(pathnames)
    file = open(pathnames)
    
    # Initialize an empty dictionary
    pathdict = Dict{Int, String}()
    
    # Read the first 4 lines of the file
    nest = 1
    for i in 1:4
        line = readline(file)
        pathdict[i] = line
    end
    commandfile = pathdict[1] * "COMMAND"
    releasefile = pathdict[1] * "RELEASES"
    outgridfile = pathdict[1] * "OUTGRID"
    input = Dict()
    command = readnamelist(commandfile)
    release = readnamelist(releasefile)
    outgrid = readnamelist(outgridfile)
    input = merge(input,command)
    input = merge(input,release)
    input = merge(input,outgrid)
    input["pathoutput"] = pathdict[2]
    input[nest] = Dict()
    input[nest]["pathmetfields"] = pathdict[3]
    input[nest]["pathavaiable"] = pathdict[4]
    readavailable!(input[nest])
    while !eof(file)
      try
        line1 = readline(file)
        line2 = readline(file)
        nest += 1
        input[nest] = Dict()
        input[nest]["pathmetfields"] = line1
        input[nest]["pathavaiable"] = line2
        readavailable!(input[nest])
      catch
        break    
      end
    end 
    input["nest"] = nest
    # Close the file
    close(file)
    
    return input
end
function readavailable!(input)
    availablefile = input["pathavaiable"]
    directory = input["pathmetfields"]
    f = open(availablefile)
    data = readlines(f)
    close(f)
    input["filename"] = []
    input["filedate"] = []
    variablesfile = split(data[1])[1]
    input["variablesfile"] = variablesfile
    for lines in data[4:end]
       entries = split(lines)
       date = DateTime(entries[1]*entries[2],dateformat"yyyymmddHHMMSS")
       filename = directory * entries[3] 
       push!(input["filename"],filename)
       push!(input["filedate"],date)
    end
    return nothing
end # end of readavailable!
function initializerun(input,metconfig)
   runconfig = Dict()
   if haskey(input,"COMMAND")
       direction = parse(Int,input["COMMAND"]["LDIRECT"][1])
       ctl = parse(Float32,input["COMMAND"]["CTL"][1])
       #get the information needed from the flexpart input files
       startdate = DateTime(input["COMMAND"]["IBDATE"][1]*input["COMMAND"]["IBTIME"][1],dateformat"yyyymmddHHMMSS")
       enddate = DateTime(input["COMMAND"]["IEDATE"][1]*input["COMMAND"]["IETIME"][1],dateformat"yyyymmddHHMMSS")
       timestep = Dates.Millisecond(input["COMMAND"]["LSYNCTIME"][1])*1000.
       outstep = Dates.Millisecond(input["COMMAND"]["LOUTSTEP"][1])*1000.
       if direction == -1 #backward
          timestep = -timestep
          outstep = -outstep
       end
       runconfig["RELEASES"] = input["RELEASE"]
       runconfig["OUTGRID"] = input["OUTGRID"]
   else
       println("need to work on the julia input data.")
       startdate = DateTime(2017,03,04,00)
       enddate =  DateTime(2017,03,04,21)
       timestep = (enddate - startdate)/8000
   end
   iout_raw = haskey(input["COMMAND"], "IOUT") ? parse(Int, input["COMMAND"]["IOUT"][1]) : 1
   runconfig["iout"] = iout_raw in (2, 10) ? 2 : iout_raw in (3, 11) ? 3 : 1
   runconfig["direction"] = direction 
   runconfig["startdate"] = startdate
   runconfig["enddate"] = enddate
   runconfig["timestep"] = timestep
   runconfig["ctl"] = ctl
   currentdate = startdate
   runconfig["currentdate"] = currentdate
   runconfig["outstep"] = outstep
   runconfig["OUTGRID"] = input["OUTGRID"]
   nest = input["nest"]
   runconfig["nest"] = nest
   runconfig["nbreleasebox"] = length(runconfig["RELEASES"]["LON1"])
   return runconfig
end
function initgridconfig(metconfig,outputconfig,nest)
   typeofgrid = metconfig["projection"]
   n = nest
   if typeofgrid == "Rectilinear"
     r_earth = 6.371e6
     local dlon::Float32 = metconfig["longitudes"][2]-metconfig["longitudes"][1]
     local dlat::Float32 = metconfig["latitudes"][2]-metconfig["latitudes"][1] 
     local dxconst::Float32 = 180. /(dlon*r_earth*pi)
     local dyconst::Float32 = 180. /(dlat*r_earth*pi)
     local Ni::Int32 = metconfig["Ni"]
     local Nj::Int32 = metconfig["Nj"]
     local Nk::Int32 = metconfig["Nk"]
     local dx::Float32 = dlon
     local dy::Float32 = dlat  
     local lonmin::Float32 = minimum(metconfig["longitudes"])
     local latmin::Float32 = minimum(metconfig["latitudes"])
     local dxout::Float32 = outputconfig["dlon"]
     local dyout::Float32 = outputconfig["dlat"]
     local xoutshift::Float32 = outputconfig[n]["xoutshift"]
     local youtshift::Float32 = outputconfig[n]["youtshift"]
     local coeffx::Float32 = dx/dxout
     local offsetx::Float32 = xoutshift/dxout
     local coeffy::Float32 = dy/dyout
     local offsety::Float32 = youtshift/dyout
     local Uedge::Bool = false
     local Vedge::Bool = false
     local Wedge::Bool = true

     gridconfig = Rectilinear( Ni, Nj, Nk, dx,dy,lonmin,latmin,dxconst, dyconst, coeffx, offsetx, coeffy, offsety, Uedge, Vedge, Wedge,dlon,dlat)
   elseif typeofgrid == "Mercator"
     r_earth_m = 6.371f6
     coslat0 = Float32(cosd(metconfig["proj_lat0"]))
     lator   = Float32(metconfig["proj_lator"])
     lonor   = Float32(metconfig["proj_lonor"])
     y_or    = Float32(r_earth_m * coslat0 * log(tand(45.0f0 + lator / 2.0f0)))
     lon0out = Float32(metconfig["longitudes"][1] - outputconfig[n]["xoutshift"])
     lat0out = Float32(metconfig["latitudes"][1]  - outputconfig[n]["youtshift"])
     gridconfig = Mercator(
        Int32(metconfig["Ni"]), Int32(metconfig["Nj"]), Int32(metconfig["Nk"]),
        Float32(metconfig["dx"]), Float32(metconfig["dy"]),
        coslat0, lonor, lator, y_or,
        lon0out, lat0out,
        Float32(outputconfig["dlon"]), Float32(outputconfig["dlat"]),
        false, false, true)
   elseif typeofgrid == "Lambert"
     local Ni2::Int32 = metconfig["Ni"]
     local Nj2::Int32 = metconfig["Nj"]
     local Nk2::Int32 = metconfig["Nk"]
     local dx2::Float32 = metconfig["dx"]
     local dy2::Float32 = metconfig["dy"]
     local dxout2::Float32 = outputconfig["dlon"]
     local dyout2::Float32 = outputconfig["dlat"]
     local xoutshift2::Float32 = outputconfig[n]["xoutshift"]
     local youtshift2::Float32 = outputconfig[n]["youtshift"]
     local proj_stdlon::Float32 = metconfig["proj_stdlon"]
     local proj_truelat1::Float32 = metconfig["proj_truelat1"]
     local proj_truelat2::Float32 = metconfig["proj_truelat2"]
     local proj_hemi::Float32 = metconfig["proj_hemi"]
     local proj_cone::Float32 = metconfig["proj_cone"]
     local proj_polei::Float32 = metconfig["proj_polei"]
     local proj_polej::Float32 = metconfig["proj_polej"]
     local proj_rebydx::Float32 = metconfig["proj_rebydx"]
     local coeffx2::Float32 = dx2/dxout2
     local offsetx2::Float32 = xoutshift2/dxout2
     local coeffy2::Float32 = dy2/dyout2
     local offsety2::Float32 = youtshift2/dyout2
     local Uedge2::Bool = false
     local Vedge2::Bool = false
     local Wedge2::Bool = true
 
     gridconfig = Lambert( Ni2, Nj2, Nk2, dx2,dy2,proj_stdlon,proj_truelat1,proj_truelat2,proj_hemi,proj_cone,proj_polei,proj_polej,proj_rebydx,coeffx2, offsetx2, coeffy2, offsety2, Uedge2, Vedge2, Wedge2)
   end
   return gridconfig
end
function initlevelconfig(metconfig,outputconfig,n)
   levelconfig = nothing
   typeoflevel = "sigmapressurelevel"
   if typeoflevel == "sigmapressurelevel"
     Nk::Int32 = metconfig["Nk"]
     ak = nothing
     bk = nothing
     try
       ak = metconfig["ak"]
       bk = metconfig["bk"] 
     catch
       ak = 1.
       bk = 0.
     end
     if metconfig["model"] == "AROME"
   akz = reverse(ak)
   bkz = reverse(bk)
   akm = zeros(Nk)
   bkm = zeros(Nk)
   akm[1] = 1.f0
   bkm[1] = 0.f0
   akm[2:end-1] = 0.5*(akz[1:end-1] + akz[2:end])
   bkm[2:end-1] = 0.5*(bkz[1:end-1] + bkz[2:end])
   akm[end] = 0.f0
   bkm[end] = 1.f0
        levelconfig = sigmapressurelevel(Nk,akm,bkm)
     else
        levelconfig = sigmapressurelevel(Nk,ak,bk)
     end
   end
   return levelconfig
end
function initarrays(metconfig,gridconfig,levelconfig)
   convertnameinput = metconfig["convertnameinput"]
   convertnameeaslee = metconfig["convertnameeaslee"]
   easlee_vars_3D = metconfig["easlee_vars_3D"]
   input_vars_3D = metconfig["input_vars_3D"]
   easlee_vars_2D = metconfig["easlee_vars_2D"]
   input_vars_2D = metconfig["input_vars_2D"]

   Ni = gridconfig.Ni
   Nj = gridconfig.Nj
   Nk = gridconfig.Nk
   if metconfig["globalx"]
     Ni += 1
   end
   metfields3D = Dict()
   metfields2D = Dict()
   metfields3Dgpu = Dict()

   if ! haskey(metconfig,"inputdimensions")
      for name in easlee_vars_3D
        metfields3D[name] = zeros(Float32,Ni,Nj,Nk,2)
      end
      for name in easlee_vars_2D
        metfields2D[name] = zeros(Float32,Ni,Nj,2)
      end
   else
      for name in easlee_vars_3D
        try  
          nameinput = convertnameeaslee[name]
          if !(nameinput == "")
            Ni,Nj,dumb = metconfig["inputdimensions"][nameinput]
            metfields3D[name] = zeros(Float32,Ni,Nj,Nk,2) #Nk has to be the staggered Nk. U,V need values at the surface
          else
            Ni,Nj,dumb = metconfig["inputdimensions"]["T"]
            metfields3D[name] = zeros(Float32,Ni,Nj,Nk,2) #Nk has to be the staggered Nk. U,V need values at the surface
          end
        catch
          continue
        end
      end
      for name in easlee_vars_2D
        try
          nameinput = convertnameeaslee[name]
          if !(nameinput == "")
           Ni,Nj = metconfig["inputdimensions"][nameinput]
           metfields2D[name] = zeros(Float32,Ni,Nj,2)
          else
            Ni,Nj = metconfig["inputdimensions"]["XLAT"]
            metfields2D[name] = zeros(Float32,Ni,Nj,2) #Nk has to be the staggered Nk. U,V need values at the surface
          end
        catch
          continue
        end
      end
   end

   # Ensure derived 3D fields always exist in the dict (they are computed by
   # metXXXfieldcalculation! but are not listed in every input variable file).
   # For Hanna/no-turbulence runs these are only accessed by advance.jl type
   # assertions; turbulence_v3(::Hanna,...) ignores them entirely, so zeros
   # are safe.  For Thomson runs they are filled correctly by the field-calc
   # routines before advance is called.
   _Ni = size(first(values(metfields3D)), 1)
   _Nj = size(first(values(metfields3D)), 2)
   _Nk = Nk
   for name in ("drhodz", "dsigwdz", "tke", "tls", "theta")
       if !haskey(metfields3D, name)
           metfields3D[name] = zeros(Float32, _Ni, _Nj, _Nk, 2)
       end
   end

   # Always pre-allocate the 2D Hanna surface fields so MetFields2D can share
   # the array reference even when the variables file doesn't list them.
   # metgribfieldcalculation! will fill them from derived quantities (stress,
   # heat flux) if they are not directly loaded from the GRIB.
   for name in ("ustar", "pblh", "wstar", "ol")
       if !haskey(metfields2D, name)
           metfields2D[name] = zeros(Float32, _Ni, _Nj, 2)
       end
   end

   return metfields3D,metfields2D,metfields3Dgpu
end


# Reset per-slot derived surface fields before loading a new met file.
# ustar/pblh/ol/wstar are DERIVED by metgribfieldcalculation! when absent
# from the file; the has_ustar/has_pblh/has_ol detection tests "slot is
# nonzero".  Without this reset, the slot retains the PREVIOUS file's derived
# values after a rotation, the detection misfires, and the surface turbulence
# parameters stay frozen at their first-load values for the whole run.
function _reset_derived_surface!(metfields2D, indt)
    # "sp" belongs here: on the ERA5 model-level path it is DERIVED from lnsp
    # only when the slot is still zero, so without this reset every pressure and
    # height in a run is built from the surface pressure of the FIRST file ever
    # loaded into the slot, and d(ps)/dt is identically zero.  Datasets that
    # archive sp reload it from the file, so zeroing first is safe there too.
    for name in ("ustar", "pblh", "wstar", "ol", "sp")
        haskey(metfields2D, name) || continue
        a = metfields2D[name]
        a isa AbstractArray || continue
        fill!(selectdim(a, ndims(a), indt), 0.0f0)
    end
    return nothing
end

# One met file → arrays slot indt: load (GRIB or NetCDF) + per-model derived
# field calculation.  Shared by readmetfiles! and the async prefetcher.
# Phase timers (decode vs derived "missing" 3D fields), for bench_grib_read.jl.
const _MET_T_DECODE = Ref(0.0)
const _MET_T_DERIVE = Ref(0.0)

function _load_met_file!(fname, indt, target_date::DateTime, metfileformat, model,
                         globalx::Bool, metconfig, metfields3D, metfields2D)
   if metfileformat == "GRIB"
      _MET_T_DECODE[] += @elapsed loadgribmetfiles!(fname, indt, false, globalx, metconfig, metfields3D, metfields2D;
                        target_date=target_date)
      # The hybrid-coordinate calculation (ak/bk hydrostatic integration) only
      # applies to ECMWF-style model-level GRIB; isobaric/height datasets have
      # their own model-specific calculation below.
      if get(metconfig, "vcoord", "hybrid") == "hybrid" && model != "AROME"
         _MET_T_DERIVE[] += @elapsed metgribfieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
      end
   elseif metfileformat == "NETCDF"
      loadnetcdfmetfiles!(fname, indt, false, target_date, metconfig, metfields3D, metfields2D)
   end
   if model == "AROME"
      metaromefieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
   elseif model == "MERRA2"
      metmerra2fieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
   elseif model == "GFS"
      metgfsfieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
   elseif model == "AROME_HDF"
      metaromehdffieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
   elseif model == "MESONH"
      metmesonhfieldcalculation!(globalx, metconfig, metfields3D, metfields2D, indt)
   else
      checkmetfields!(metconfig, metfields3D, metfields2D, indt)
   end
   return nothing
end

"""
    readmetfiles!(nestctx, currentdate, direction, metconfig,
                  metfields3D, metfields2D, metfields3Dgpu)

Keep the two met time slots bracketing `currentdate`: on the first call load
both files, afterwards rotate in the next file (in the run `direction`) when
the current date passes the newer slot.  FLEXPART lineage: getfields/readwind.
"""
function readmetfiles!(nestctx, currentdate, direction, metconfig, metfields3D, metfields2D, metfields3Dgpu)
   # ── Typed inputs ──────────────────────────────────────────────────────────────
   n           = nestctx.id
   mc          = nestctx.metconfig          # MetConfig typed struct
   filedate    = mc.filedate
   filename    = mc.filename
   metfileformat = mc.metfileformat
   model       = mc.model
   globalx     = nestctx.gridconfig.globalx

   # Mutable local t_index vector — sub-functions mutate it in-place, then we
   # write the final bracket back to the typed NestRunConfig at the end.
   t_idx = Int[nestctx.nestrunconfig.t_index[1], nestctx.nestrunconfig.t_index[2]]

   #load the metfields
   ###############################
   first_call = (t_idx[1] == 0 && t_idx[2] == 0)
   if first_call  #need to read 2 files for the first time
      # review 2026-09-04 A3 (item 3): bracket clamped to [1, nfile] and an explicit
      # error when the start lies outside the archive (was BoundsError / MethodError).
      t_index = _first_met_bracket(filedate, currentdate, direction)
      t_idx[1] = t_index
      fname = filename[t_index]
      _reset_derived_surface!(metfields2D, 1)
      println("index=",1,"file to read=",fname)
      _load_met_file!(fname, 1, filedate[t_idx[1]], metfileformat, model,
                      globalx, metconfig, metfields3D, metfields2D)
      t_index = direction == 1 ? t_index + 1 : t_index - 1
      t_idx[2] = t_index
      fname = filename[t_index]
      _reset_derived_surface!(metfields2D, 2)
      println("index=",2,"file to read=",fname)
      _load_met_file!(fname, 2, filedate[t_idx[2]], metfileformat, model,
                      globalx, metconfig, metfields3D, metfields2D)
   else #need to read a new file?
      if direction == 1
         t_index = maximum(t_idx)
         testdate = currentdate > filedate[t_index]
      else
         t_index = minimum(t_idx)
         testdate = currentdate < filedate[t_index]
      end
      # Would the next bracket step read PAST the archive?  Forward: max+1 must
      # exist; backward: min−1 must exist.  If not, hold the current (boundary)
      # bracket — the run continues on the last available met (met_time_index_array
      # clamps the weight to the frame), instead of an out-of-bounds file read.
      nfile = length(filedate)
      next_ok = direction == 1 ? maximum(t_idx) + 1 <= nfile :
                                 minimum(t_idx) - 1 >= 1
      if testdate && !next_ok
         println("nothing to read (end of met archive — holding the last frame)")
         testdate = false
      end
      if testdate
         fname = nothing
         indt  = nothing
         if t_idx[1] < t_idx[2]
            if direction == 1
               t_idx[1] = t_idx[2] + 1   # new data go into slot 1
               indt  = 1
               fname = filename[t_idx[1]]
            else
               t_idx[2] = t_idx[1] - 1   # new data go into slot 2
               indt  = 2
               fname = filename[t_idx[2]]
            end
         else
            if direction == 1
               t_idx[2] = t_idx[1] + 1   # new data go into slot 2
               indt  = 2
               fname = filename[t_idx[2]]
            else
               t_idx[1] = t_idx[2] - 1   # new data go into slot 1
               indt  = 1
               fname = filename[t_idx[1]]
            end
         end
         _reset_derived_surface!(metfields2D, indt)
         println("index=",indt,"file to read=",fname)
         _load_met_file!(fname, indt, filedate[t_idx[indt]], metfileformat, model,
                         globalx, metconfig, metfields3D, metfields2D)
      else
         println("nothing to read")
      end
   end
   # Write the updated time-bracket indices back to the typed NestRunConfig.
   nestctx.nestrunconfig.t_index = (t_idx[1], t_idx[2])
   return nothing
end
# review 2026-09-04 A3 (item 3): index of the FIRST file of the initial two-file
# bracket (the second is +1 forward, −1 backward), clamped so both indices are
# valid.  Errors, naming the run start and the archive span, when the start lies
# outside [filedate[1], filedate[end]] — the legacy loop gave t_index = 0/nfile+1
# (BoundsError) or `nothing` (MethodError) there.  filedate must be ascending.
function _first_met_bracket(filedate::AbstractVector, currentdate::DateTime, direction::Integer)
   nfile = length(filedate)
   nfile >= 2 || error("met archive holds $(nfile) file(s); at least two analysis times are needed")
   (filedate[1] <= currentdate <= filedate[end]) ||
      error("run start $(currentdate) lies outside the met archive span " *
            "[$(filedate[1]) … $(filedate[end])] ($(nfile) files, direction=$(direction))")
   if direction == 1
      i = findfirst(>=(currentdate), filedate)
      return clamp(i - 1, 1, nfile - 1)
   else
      i = findlast(<=(currentdate), filedate)
      return clamp(i + 1, 2, nfile)
   end
end

# review 2026-09-04 A3 (item 6): per-file bookkeeping of the GRIB message loop —
# levels stored per 3-D input variable, plus the accumulation step of the sshf
# message (for item 1).  Atomics: the reader runs several decode workers.
struct GribScan
   hits::Dict{String,Threads.Atomic{Int}}   # 3-D input shortName → levels stored
   flux_step_s::Threads.Atomic{Int}         # endStep−startStep of sshf [s]; −1 = unseen
end
GribScan(input_vars_3D) = GribScan(
   Dict{String,Threads.Atomic{Int}}(String(n) => Threads.Atomic{Int}(0) for n in input_vars_3D),
   Threads.Atomic{Int}(-1))

# A message set that is incomplete (companion file absent, truncated file, level
# subset) leaves the PREVIOUS file's data in the rotated slot silently.  The
# first file of the run sets the per-variable reference level count; later files
# are compared against it.  Zero levels for T/U/V is fatal.
function _check_grib_completeness!(metconfig, scan::GribScan, filename, convertname, levelmap, Nk)
   expected = levelmap === nothing ? Int(Nk) - 1 : length(levelmap)
   first = !haskey(metconfig, "grib_nlev3D")
   ref   = get!(metconfig, "grib_nlev3D", Dict{String,Int}())
   for (name, a) in scan.hits
      n = a[]
      ename = string(get(convertname, name, name))
      if n == 0 && ename in ("T", "Ux", "Vy")
         error("GRIB read: no '$name' ($ename) message in $(basename(filename)) — required " *
               "3-D field missing (truncated file, wrong companion pattern or level type)")
      end
      if first
         n > 0 && (ref[name] = n)
         0 < n < expected && @warn "GRIB read: 3-D field has fewer levels than the ladder" file=basename(filename) var=name easlee=ename levels_found=n levels_expected=expected
      elseif haskey(ref, name) && n != ref[name]
         @warn "GRIB read: incomplete 3-D field — slot keeps the PREVIOUS file's levels" file=basename(filename) var=name easlee=ename levels_found=n levels_expected=ref[name]
      end
   end
   # accumulation period of the surface fluxes (item 1), recorded once per run
   if !haskey(metconfig, "sfcflux_step_s") && scan.flux_step_s[] >= 0
      metconfig["sfcflux_step_s"] = scan.flux_step_s[]
   end
   return nothing
end

function loadgribmetfiles!(filename,indt,flipit,globalx::Bool,metconfig,metfields3D,metfields2D;
                           target_date::Union{Nothing,DateTime}=nothing)
   convertname = metconfig["convertnameinput"]
   easlee_vars_3D = metconfig["easlee_vars_3D"]
   input_vars_3D = metconfig["input_vars_3D"]
   easlee_vars_2D = metconfig["easlee_vars_2D"]
   input_vars_2D = metconfig["input_vars_2D"]

   Nk = metconfig["Nk"]
   Ni = metconfig["Ni"]
   # Level mapping for non-hybrid GRIB datasets (isobaric / height levels):
   # message level VALUE → bottom-up ladder index built at gridcheck time.
   levelmap = nothing
   leveltol = ""
   if get(metconfig, "vcoord", "hybrid") != "hybrid"
      levelmap = Dict{Float64,Int}(v => k for (k, v) in enumerate(metconfig["levelvalues"]))
      leveltol = metconfig["leveltol"]
   end
   # Multi-step files (several valid times per GRIB): filter by validity time.
   tdate = get(metconfig, "multistep", false) ? target_date : nothing
   crop  = get(metconfig, "crop", nothing)     # domain crop window (_crop_metconfig!)
   flip_lat = get(metconfig, "flip_lat", true)::Bool   # review 2026-09-04 A3 (item 5)
   scan  = GribScan(input_vars_3D)                      # review 2026-09-04 A3 (item 6)
   _read_grib_messages!(filename, indt, flipit, globalx, Ni, Nk,
                        convertname, input_vars_3D, input_vars_2D,
                        metfields3D, metfields2D;
                        levelmap=levelmap, leveltol=leveltol, target_date=tdate,
                        crop=crop, flip_lat=flip_lat, scan=scan)
   # Companion files (datasets that split fields across several files at one
   # valid time, e.g. ERA5 IFS_EA_MET_/IFS_EA_CLD_ or AROME-HDF HP1/SP1).
   # Patterns come from the namelist "#info companions" entry; only variables
   # listed in the variables file are picked up.
   for (pfrom, pto) in get(metconfig, "companions", Tuple{String,String}[])
      occursin(pfrom, basename(filename)) || continue
      cname = joinpath(dirname(filename), replace(basename(filename), pfrom => pto))
      println("  companion: ", basename(cname), isfile(cname) ? "" : "  (MISSING)")
      if isfile(cname)
         _read_grib_messages!(cname, indt, flipit, globalx, Ni, Nk,
                              convertname, input_vars_3D, input_vars_2D,
                              metfields3D, metfields2D;
                              levelmap=levelmap, leveltol=leveltol, target_date=tdate,
                              crop=crop, flip_lat=flip_lat, scan=scan)
      end
   end
   _check_grib_completeness!(metconfig, scan, filename, convertname, levelmap, Nk)
   # Pressure-level stitch (AROME-HDF IP1): read the isobaric package's u/v/t/r/z
   # onto the UPPER ladder slots n_hp+2 … Nk via an isobaric level map.  These
   # slots sit above the height block and are left untouched by the surface shift
   # in metaromehdffieldcalculation!, so they land directly at their final index.
   pc = get(metconfig, "pressure_companion", nothing)
   if pc !== nothing && haskey(metconfig, "ip_pressures") && occursin(pc[1], basename(filename))
      ipfile = joinpath(dirname(filename), replace(basename(filename), pc[1] => pc[2]))
      if isfile(ipfile)
         n_hp   = metconfig["n_hp"]::Int
         ippres = metconfig["ip_pressures"]::Vector{Float64}
         ip_levelmap = Dict{Float64,Int}(p => n_hp + 1 + r for (r, p) in enumerate(ippres))
         println("  pressure-stitch: ", basename(ipfile))
         _read_grib_messages!(ipfile, indt, flipit, globalx, Ni, Nk,
                              convertname, input_vars_3D, input_vars_2D,
                              metfields3D, metfields2D;
                              levelmap=ip_levelmap, leveltol="isobaricInhPa", target_date=tdate,
                              crop=crop, flip_lat=flip_lat)
      end
   end
   return nothing
end #loadmetconfig

# Decode one GRIB message (already wrapped) and scatter it into the met arrays.
# Pure per-message work: each message targets a distinct (field, level, indt) slice,
# so concurrent calls on different messages never write the same cell — safe to run
# from several tasks (the `values` unpack is the costly part we want parallel).
# review 2026-09-04 A3 (item 7): decode → Float32 slot in ONE pass — the
# missing-value test inline and the latitude flip as an index reversal — instead
# of two full Float64 temporaries per message (`map(x -> x == mv ? NaN : x, v)`
# then `reverse(v, dims=2)`).  Bit-identical to the old store (same Float64→
# Float32 rounding, NaN → NaN32).  `dst` is the (Ni, Nj) slot view; `crop` is
# the (iidx, jr) gather window of _crop_metconfig! or nothing.
function _store_grib_slice!(dst::AbstractMatrix{Float32}, vals::AbstractMatrix{Float64},
                            has_mv::Bool, mv::Float64, flip_lat::Bool, crop)
   Nj = size(vals, 2)
   if crop === nothing
      Ni = size(dst, 1)
      size(dst, 2) == Nj || error("GRIB message has $(Nj) rows, met slot has $(size(dst, 2))")
      @inbounds for j in 1:Nj
         js = flip_lat ? Nj - j + 1 : j
         for i in 1:Ni
            x = vals[i, js]
            dst[i, j] = (has_mv && x == mv) ? NaN32 : Float32(x)
         end
      end
   else
      iidx, jr = crop
      @inbounds for (j, jd) in enumerate(jr)
         js = flip_lat ? Nj - jd + 1 : jd
         for (i, id) in enumerate(iidx)
            x = vals[id, js]
            dst[i, j] = (has_mv && x == mv) ? NaN32 : Float32(x)
         end
      end
   end
   return nothing
end

function _apply_grib_message!(test, indt, flipit, globalx::Bool, Ni, Nk,
                              convertname, input_vars_3D, input_vars_2D,
                              metfields3D, metfields2D, levelmap, leveltol, target_date,
                              crop = nothing, flip_lat::Bool = true,
                              scan::Union{Nothing,GribScan} = nothing)
   name = test["shortName"]
   is3D = name in input_vars_3D
   is2D = name in input_vars_2D
   (is3D || is2D) || return nothing
   if target_date !== nothing
      vt = Int(test["validityTime"])
      msgdate = DateTime(string(test["validityDate"]), dateformat"yyyymmdd") +
                Dates.Hour(div(vt, 100)) + Dates.Minute(mod(vt, 100))
      msgdate == target_date || return nothing
   end
   # Mask bitmap-missing values (eccodes returns 9999 by default, e.g. AROME-HDF
   # corners outside the native Lambert domain).  NaN is handled by the physics kernels.
   msgvalues = test["values"]::Matrix{Float64}
   has_mv = Bool(test["bitmapPresent"])
   mv = has_mv ? Float64(test["missingValue"]) : NaN
   if is3D
      local level::Int
      if levelmap === nothing
         level = Nk - test["level"]
      else
         # Non-hybrid: keep only the namelist's level type (e.g. GFS also carries
         # isobaricInPa strato duplicates) and map the level value onto the ladder index.
         test["typeOfLevel"] == leveltol || @goto twod
         level = get(levelmap, Float64(test["level"]), 0)
         level == 0 && @goto twod
      end
      if flipit
         metfields3D[convertname[name]][:,:,level,1] = metfields3D[convertname[name]][:,:,level,2]
      end
      # domain crop: store only the box+halo sub-window of the decoded field
      # (iidx may roll a seam-crossing window contiguous; jr slices latitude)
      _store_grib_slice!(view(metfields3D[convertname[name]], 1:Ni, :, level, indt),
                         msgvalues, has_mv, mv, flip_lat, crop)
      scan === nothing || Threads.atomic_add!(scan.hits[name], 1)   # review 2026-09-04 A3 (item 6)
      if globalx
         metfields3D[convertname[name]][Ni+1,:,level,indt] = metfields3D[convertname[name]][1,:,level,indt]
      end
   end
   @label twod
   if is2D
      # A 2D var must not swallow 3D messages sharing its shortName — but only when
      # the message level is actually on the 3D ladder (e.g. AROME-HDF 2t/2r at 2 m
      # are heightAboveGround like the 3D fields, yet legitimate 2D variables).
      (levelmap !== nothing && test["typeOfLevel"] == leveltol &&
       haskey(levelmap, Float64(test["level"]))) && return nothing
      if flipit
         metfields2D[convertname[name]][:,:,1] = metfields2D[convertname[name]][:,:,2]
      end
      _store_grib_slice!(view(metfields2D[convertname[name]], 1:Ni, :, indt),
                         msgvalues, has_mv, mv, flip_lat, crop)
      # review 2026-09-04 A3 (item 1): accumulation step of the heat-flux message
      # (endStep − startStep, in seconds); 0 on de-accumulated flex_extract files.
      if scan !== nothing && convertname[name] == "sshf"
         try
            su = Int(test["stepUnits"])                    # 0 = min, 1 = h, 2 = day
            fac = su == 0 ? 60 : su == 2 ? 86400 : 3600
            scan.flux_step_s[] = (Int(test["endStep"]) - Int(test["startStep"])) * fac
         catch
         end
      end
      if globalx
         metfields2D[convertname[name]][Ni+1,:,indt] = metfields2D[convertname[name]][1,:,indt]
      end
   end
   return nothing
end

# GRIB read.  Decoding (unpacking `values`) dominates for packed ERA5/IFS fields, so it
# is parallelised the FLEXPART way: each worker opens its OWN file handle and decodes a
# strided subset of the messages (worker t handles messages t, t+nt, …).  The values
# unpack runs concurrently; each message targets a disjoint (field,level,indt) slice so
# the writes never collide.  GRIB_NT caps the worker count (default 4 — decode is memory-
# bandwidth bound, so more rarely helps and competes with the rest of the run); GRIB_NT=1
# keeps the original single-handle serial reader.
function _read_grib_messages!(filename, indt, flipit, globalx::Bool, Ni, Nk,
                              convertname, input_vars_3D, input_vars_2D,
                              metfields3D, metfields2D;
                              levelmap=nothing, leveltol::String="",
                              target_date::Union{Nothing,DateTime}=nothing,
                              crop=nothing, flip_lat::Bool=true,
                              scan::Union{Nothing,GribScan}=nothing)
   nt = clamp(parse(Int, get(ENV, "GRIB_NT", "4")), 1, max(1, Threads.nthreads()))

   if nt <= 1
      GribFile(filename) do f
         for msg in f
            _apply_grib_message!(msg, indt, flipit, globalx, Ni, Nk,
                                 convertname, input_vars_3D, input_vars_2D,
                                 metfields3D, metfields2D, levelmap, leveltol, target_date,
                                 crop, flip_lat, scan)
         end
      end
      return nothing
   end

   # Warm eccodes' lazy definition cache single-threaded (one message) before the
   # workers touch the shared default context concurrently.
   GribFile(filename) do f
      for msg in f
         msg["shortName"]; try; msg["values"]; catch; end
         break
      end
   end

   @sync for t in 1:nt
      Threads.@spawn GribFile(filename) do f
         c = 0
         for msg in f
            c += 1
            ((c - 1) % nt == (t - 1)) || continue
            _apply_grib_message!(msg, indt, flipit, globalx, Ni, Nk,
                                 convertname, input_vars_3D, input_vars_2D,
                                 metfields3D, metfields2D, levelmap, leveltol, target_date,
                                 crop, flip_lat, scan)
         end
      end
   end
   return nothing
end
# Staging buffer for in-place NetCDF variable reads (loadnetcdfmetfiles!):
# one Float32 buffer per variable SHAPE, cached in metconfig and reused across
# variables and files — replaces the per-read `v.vars[...]` + `reverse(...)`
# temporaries (two full-array allocations per variable per file).
# Always returns a 4-D array (Ni,Nj,Nk,1); 3-D file variables are read through
# a reshaped view of the same memory.
function _ncf_stage!(metconfig, name::AbstractString, dims::Tuple)
    cache = get!(metconfig, "_ncstage", Dict{Tuple,Array{Float32,4}}())
    d3 = length(dims) >= 3 ? (Int(dims[1]), Int(dims[2]), Int(dims[3])) :
                             (Int(dims[1]), Int(dims[2]), 1)
    get!(cache, d3) do
        Array{Float32,4}(undef, d3[1], d3[2], d3[3], 1)
    end
end

function loadnetcdfmetfiles!(filename,indt,flipit,target_date::DateTime,metconfig,metfields3D,metfields2D)
   _read_netcdf_file!(filename, indt, target_date, metconfig, metfields3D, metfields2D)
   # Companion files (e.g. MERRA2 asm_Nv → mst_Ne / trb_Ne collections); see
   # the namelist "#info companions" entry.
   for (pfrom, pto) in get(metconfig, "companions", Tuple{String,String}[])
      occursin(pfrom, basename(filename)) || continue
      cname = joinpath(dirname(filename), replace(basename(filename), pfrom => pto))
      isfile(cname) && _read_netcdf_file!(cname, indt, target_date, metconfig, metfields3D, metfields2D)
   end
   return nothing
end
# review 2026-09-04 A3 (item 4): index of the WRF `Times` slice equal to `filedate`.
function _wrf_time_slice(Times::AbstractMatrix, filedate::DateTime, filename)
   dates = [DateTime(join(view(Times, :, k)), dateformat"yyyy-mm-dd_HH:MM:SS") for k in 1:size(Times, 2)]
   k = findfirst(==(filedate), dates)
   k === nothing && error("WRF file $(basename(string(filename))) has no time slice for $(filedate); " *
                          "the file holds $(length(dates)) slice(s) spanning $(dates[1]) … $(dates[end])")
   return k
end
function _read_netcdf_file!(filename, indt, filedate::DateTime, metconfig, metfields3D, metfields2D)
   convertname = metconfig["convertnameinput"]
   input_vars_3D = metconfig["input_vars_3D"]
   input_vars_2D = metconfig["input_vars_2D"]
   Nk = metconfig["Nk"]
   model = metconfig["model"]
   v = NetCDF.open(filename)
   #it is assumed here that a single filename can have multiple dates.
   inct = 1
   if model == "WRF"
      # review 2026-09-04 A3 (item 4): search the TIME dimension (columns of the
      # 19×ntime char array), not the 19 characters of the first date string, and
      # error when no slice matches instead of silently using the last one.
      inct = _wrf_time_slice(v.vars["Times"][:,:], filedate, filename)
   elseif model == "MERRA2"
      # time(time) in "minutes since YYYY-MM-DD HH:MM:SS"; pick nearest entry.
      tvar = v.vars["time"]
      m = match(r"since (\d{4}-\d{2}-\d{2})[ T](\d{2}:\d{2}:\d{2})", tvar.atts["units"])
      tbase = DateTime(m.captures[1] * "T" * m.captures[2])
      tvals = tvar[:]
      inct = argmin([abs(Dates.value(tbase + Dates.Minute(round(Int, t)) - filedate)) for t in tvals])
   end
   if metconfig["model"] == "AROME"
      inct = 1
   end
   has_dims = haskey(metconfig, "inputdimensions")
   for name in input_vars_3D
      haskey(v.vars, name) || continue           # var lives in another companion file
      if has_dims
      dumb,dumb,Nk2 = metconfig["inputdimensions"][name]
      # Reuse ONE staging buffer per variable shape (was: v.vars[...] alloc + reverse
      # alloc = two full temporaries per var per file → GiB-scale GC churn, the top
      # cost of the met read).  NetCDF reads in-place into the staging buffer; the
      # k-flip (AROME stores top-down) is a direct indexed copy into the slot.
      stage = _ncf_stage!(metconfig, name, size(v.vars[name]))
      if ndims(v.vars[name]) == 4
         NetCDF.readvar!(v, name, stage; start=[1,1,1,inct], count=[size(stage,1),size(stage,2),size(stage,3),1])
      else
         NetCDF.readvar!(v, name, stage; start=[1,1,1], count=[size(stage,1),size(stage,2),size(stage,3)])
      end
      kupper = Nk2 < Nk ? Nk-1 : Int(Nk)
      dst3 = metfields3D[convertname[name]]
      nk_src = size(stage, 3)
      if metconfig["model"] == "AROME"
         @inbounds for k in 1:min(kupper, nk_src)
            ksrc = nk_src - k + 1                      # reverse(dims=3)
            @views dst3[:, :, k, indt] .= stage[:, :, ksrc, 1]
         end
      else
         @inbounds for k in 1:min(kupper, nk_src)
            @views dst3[:, :, k, indt] .= stage[:, :, k, 1]
         end
      end
      else
         # Generic path (MERRA2, MESONH): uniform grid, allocate-by-Nk arrays.
         var  = v.vars[name]
         data = ndims(var) == 4 ? var[:,:,:,inct] : var[:,:,:]
         Nk2  = size(data, 3)
         Ni0  = size(data, 1)
         dst  = metfields3D[convertname[name]]
         if model == "MERRA2"
            # MERRA2 levels are top-down (lev 1 = model top): flip to the
            # EASLEE bottom-up convention.  The +1 surface-slot shift is done
            # in metmerra2fieldcalculation!.
            dst[1:Ni0,:,1:Nk2,indt] = reverse(Float32.(data), dims=3)
         else
            dst[1:Ni0,:,1:min(Nk2,Int(Nk)),indt] = Float32.(data[:,:,1:min(Nk2,Int(Nk))])
         end
         # Global longitude wrap column (arrays carry Ni+1 points).
         size(dst, 1) == Ni0 + 1 && (dst[Ni0+1,:,:,indt] = dst[1,:,:,indt])
      end
   end
   for name in input_vars_2D
      haskey(v.vars, name) || continue
      var = v.vars[name]
      dst = metfields2D[convertname[name]]
      Ni0 = size(var, 1)
      if ndims(var) == 3
         dst[1:Ni0,:,indt] = Float32.(var[:,:,inct])
      else
         # Static 2D field (orography, land mask, MESONH ZS/LAT/LON…)
         dst[1:Ni0,:,indt] = Float32.(var[:,:])
      end
      size(dst, 1) == Ni0 + 1 && (dst[Ni0+1,:,indt] = dst[1,:,indt])
   end
end
function checkmetfields!(metconfig,metfields3D,metfields2D,indt)
# Check that the derived 3D fields (zlev, Ux, Vy, levdot, rho) are defined;
# if not, try to calculate them from the fields that have been loaded
# (generic NetCDF path — mainly WRF-style inputs).
   convertnameeaslee = metconfig["convertnameeaslee"]
   convertname = metconfig["convertnameinput"]
   easlee_vars_3D = metconfig["easlee_vars_3D"]
   input_vars_3D = metconfig["input_vars_3D"]
   easlee_vars_2D = metconfig["easlee_vars_2D"]
   input_vars_2D = metconfig["input_vars_2D"]
#3D:   zlev
   if convertnameeaslee["zlev"] == ""
      #need to calculate height from geopotential
      ph = 0.
      if haskey(metfields3D,"ph")
         phb = 0.
         ph = metfields3D["ph"][:,:,:,indt]
         if haskey(metfields3D,"phb")
            phb = metfields3D["phb"][:,:,:,indt]
         end
         ph += phb
      # (deriving heights from ak/bk instead, when ph is absent, is not implemented)
      metfields3D["zlev"][:,:,:,indt] = ph/9.81
      surfacez = metfields3D["zlev"][:,:,1,indt]
      metfields3D["zlev"][:,:,:,indt] = metfields3D["zlev"][:,:,:,indt] .- reshape(surfacez,size(surfacez, 1), size(surfacez, 2), 1)
      end
   end
#3D:   Ux. Horizontal wind along x axis 
   if convertnameeaslee["Ux"] == ""
     stoprun = true
     #such a case cannot happen on a regular lat/lon grid, as zonal wind and Ux are equal
     if metconfig["projection"] == "Rectilinear"
        print("The model has to stop, as zonal wind is missing")
     else
        if "AVGFLX_RUM" in metconfig["easlee_vars_3D"]
           #we need to fix the time of the met field, because it is shifted by the WRF output time step/2.
           #for the moment, we will assume it is ok
           mu = metfields2D["muu"][:,:,indt]
           # AVGFLX_RUM is on the U-stagger grid (Ni+1 points); destagger to mass
           # points (Ni) by averaging adjacent stagger-point values before dividing
           # by mu (which lives on mass points).
           rum_mass = 0.5f0 .* (metfields3D["AVGFLX_RUM"][1:end-1,:,:,indt] .+
                                 metfields3D["AVGFLX_RUM"][2:end,  :,:,indt])
           metfields3D["Ux"][:,:,:,indt] = rum_mass ./ reshape(mu, size(mu, 1), size(mu, 2), 1)
           stoprun = false
        elseif "Uzonal" in metconfig["easlee_vars_3D"]
           #need to correct it along the x axis using map factors
           if "MAPFAC_UX" in metconfig["easlee_vars_2D"]
              m_u = metfields2D["MAPFAC_UX"][:,:,1,indt]
              metfields3D["Ux"][:,:,indt] = metfields3D["Uzonal"][:,:,indt]./m_u # need to fix this calculation
              stoprun = false
           elseif "MAPFAC_MX" in metconfig["easlee_vars_2D"]
              m_x = metfields2D["MAPFAC_MX"][:,:,1,indt]
      if size(m_x[:,:,1,indt]) == size(metfields3D["Uzonal"][:,:,1,1])
	#check that it has the same size as Uzonal or not. If yes, then Uzonal is on the center grid
                metfields3D["Ux"][:,:,:,indt] = metfields3D["Uzonal"][:,:,:,indt] ./ m_x # need to fix this calculation
                stoprun = false
              else
                #Uzonal is on the edges and MAPFAC is at the center grid. No choice but to interpolate MAPFAC
                print("MAPFAC is at the center of grids while U is on the edges. Need to fix the problem")
                stoprun = true #will need to find a fix later
              end
           end
        end
     end
   end
     #assume that we have zonal wind on a lambert conformal grid
#3D:   Vy
   if convertnameeaslee["Vy"] == ""
     stoprun = true
     #such a case cannot happen on a regular lat/lon grid, as meridional wind and Vy are equal
     if metconfig["projection"] == "Rectilinear"
        print("The model has to stop, as meridional wind is missing")
     else
        if "AVGFLX_RVM" in metconfig["easlee_vars_3D"]
           #we need to fix the time of the met field, because it is shifted by the WRF output time step/2.
           #for the moment, we will assume it is ok
           mu = metfields2D["muu"][:,:,indt]   # total dry air column mass (same as Ux)
           # AVGFLX_RVM is on the V-stagger grid (Nj+1 points); destagger to mass
           # points (Nj) by averaging adjacent stagger-point values.
           rvm_mass = 0.5f0 .* (metfields3D["AVGFLX_RVM"][:,1:end-1,:,indt] .+
                                 metfields3D["AVGFLX_RVM"][:,2:end,  :,indt])
           metfields3D["Vy"][:,:,:,indt] = rvm_mass ./ reshape(mu, size(mu, 1), size(mu, 2), 1)
           stoprun = false
        elseif "Vmerid" in metconfig["easlee_vars_3D"]
           #need to correct it along the x axis using map factors
           if "MAPFAC_Vy" in metconfig["easlee_vars_2D"]
              m_v = metfields2D["MAPFAC_VY"][:,:,1,indt]
              metfields3D["Vy"][:,:,:,indt] = metfields3D["Vmerid"][:,:,:,indt]./m_v # need to fix this calculation
              stoprun = false
           elseif "MAPFAC_My" in metconfig["easlee_vars_2D"]
              m_y = metfields2D["MAPFAC_MY"][:,:,1,indt]
              if size(m_y[:,:,1]) == size(metfields3D["Vmerid"][:,:,1,1])
                #check that it has the same size as Vmerid or not. If yes, then Vmerid is on the center grid
                metfields3D["Vy"][:,:,:,indt] = metfields3D["Vmerid"][:,:,:,indt]/m_v # need to fix this calculation
                stoprun = false
              else
                #Vmerid is on the edges and MAPFAC is at the center grid. No choice but to interpolate MAPFAC
                print("MAPFAC is at the center of grids while U is on the edges. Need to fix the problem")
                stoprun = true #will need to find a fix later
              end
           end
        end
     end
   end

#3D:   levdot
   if convertnameeaslee["levdot"] == ""
     stoprun = true
     if metconfig["projection"] == "Rectilinear"
        print("The model has to stop, as meridional wind is missing")
     else
        if "AVGFLX_WWM" in metconfig["easlee_vars_3D"]
           #we need to fix the time of the met field, because it is shifted by the WRF output time step/2.
           #for the moment, we will assume it is ok
           mu   = metfields2D["muu"][:,:,indt]   # total dry air column mass (same as Ux/Vy)
           # levdot currently holds raw AVGFLX_WWM on w-stagger (Nk_w = Nk_mass+1 levels).
           # Infer mass-level count from the loaded array rather than metconfig["Nk"],
           # which counts w-stagger levels for WRF.
           Nk_w    = size(metfields3D["levdot"], 3)   # w-stagger level count (e.g. 51)
           Nk_mass = Nk_w - 1                          # mass level count (e.g. 50)
           # Destagger from w-stagger to mass levels, then convert etadot [1/s] to
           # dindl/dt [level/s]: multiply by (Nk_mass-1) and negate (WRF eta increases
           # downward while EASLEE indl increases upward).
           wwm_mass = 0.5f0 .* (metfields3D["levdot"][:,:,1:end-1,indt] .+
                                 metfields3D["levdot"][:,:,2:end,  indt])
           metfields3D["levdot"][:,:,1:Nk_mass,indt] = -wwm_mass ./ reshape(mu, size(mu,1), size(mu,2), 1) .* Float32(Nk_mass - 1)
           stoprun = false
        elseif "W" in metconfig["easlee_vars_3D"] 
            zlev = metfields3D["zlev"][:,:,:,indt]
            Nk = metconfig["Nk"]
            metfields3D["levdot"][:,:,2:Nk-1,indt] = metfields3D["W"][:,:,2:Nk-1,indt] ./ (zlev[:,:,3:Nk]-zlev[:,:,1:Nk-2])/2.f0
            metfields3D["levdot"][:,:,1,indt] = metfields3D["W"][:,:,1,indt] ./ (zlev[:,:,2]-zlev[:,:,1])
            metfields3D["levdot"][:,:,Nk,indt] = metfields3D["W"][:,:,Nk,indt] ./ (zlev[:,:,Nk]-zlev[:,:,Nk-1])
        end
     end
   end 

#3D:   rho
   if convertnameeaslee["rho"] == ""
      #need to calculate rho from T and P
      P = 0.
      if haskey(metfields3D,"p")
         PB = 0.
         P = metfields3D["p"][:,:,:,indt]
         if haskey(metfields3D,"pb")
            PB = metfields3D["pb"][:,:,:,indt]
         end
         P += PB
      # (deriving P from ak/bk instead, when p is absent, is not implemented)
      end
      if metconfig["model"] == "WRF"
         # WRF stores T as perturbation potential temperature θ′ = θ − 300 K.
         # Recover actual temperature via the Exner function:
         #   θ = θ′ + 300,   T = θ × (P/P₀)^(R_d/c_p) with P₀=100000 Pa, R_d/c_p=0.2854
         theta = metfields3D["T"][:,:,:,indt] .+ 300.0f0
         T = theta .* (P ./ 100000.0f0) .^ 0.2854f0
      else
         T = metfields3D["T"][:,:,:,indt]
      end
      metfields3D["rho"][:,:,:,indt] = P ./ T ./ 287.05
   end
   
end
# ── Hanna (1982) / FLEXPART helper functions ──────────────────────────────────
# Duplicated here so initialize.jl (Main scope) can call hanna_fill_tke!
# without depending on the Operators module.  The authoritative implementation
# lives in src/physics/turbulence_hanna.jl.
@inline function _hanna_wstar(ustar::Float32, pblh::Float32, ol::Float32)
    ol < 0.0f0 || return 0.0f0
    return ustar * cbrt(max(-pblh / ol, 0.0f0))
end

# Returns (sigw [m/s], tlw [s]) using FLEXPART stability regimes, with a
# blended entrainment zone [0.85h, 1.15h] joining the BL and free-troposphere
# values (no hard switch at h — see hanna_sigw_tlw in turbulence_hanna.jl,
# the authoritative copy, for the rationale).
@inline function _hanna_sigw_tlw(z::Float32, pblh::Float32, ustar::Float32,
                                   wstar::Float32, ol::Float32, sigw_ft::Float32)
    h  = max(pblh, 50.0f0)
    z1 = 0.85f0 * h
    z2 = 1.15f0 * h
    sigw_f = max(sigw_ft, 0.05f0)
    tlw_f  = clamp(100.0f0 * sigw_f, 0.5f0, 300.0f0)
    z >= z2 && return sigw_f, tlw_f
    sigw_b, tlw_b = _hanna_bl_main(z, pblh, ustar, wstar, ol)
    z <= z1 && return sigw_b, tlw_b
    x = (z - z1) / (z2 - z1)
    w = 0.5f0 * (1.0f0 + cos(Float32(π) * x))
    sigw = sqrt(max(w * sigw_b^2 + (1.0f0 - w) * sigw_f^2, 0.25f0 * sigw_f^2))
    tlw  = clamp(w * tlw_b + (1.0f0 - w) * tlw_f, 0.5f0, 300.0f0)
    return sigw, tlw
end

# Hanna BL regime formulas (sigw, tlw only — drift gradients are computed
# online by the Thomson kernel from the filled tke profile).
@inline function _hanna_bl_main(z::Float32, pblh::Float32, ustar::Float32,
                                 wstar::Float32, ol::Float32)
    h   = max(pblh, 50.0f0)
    aol = abs(ol)
    # Neutral: pblh/|L| < 1
    if aol < 1.0f-6 || h / aol < 1.0f0
        ust   = max(ustar, 1.0f-4)
        corr  = z / ust
        sigw0 = 1.3f0 * ust * exp(-2.0f-4 * corr)
        sigw  = sigw0 + 1.0f-2
        tlw   = clamp(0.5f0 * z / max(sigw, 1.0f-6) /
                      max(1.0f0 + 1.5f-3 * corr, 1.0f-6), 0.5f0, 300.0f0)
    # Convective: L < 0, pblh/|L| ≥ 1
    elseif ol < 0.0f0
        zeta = max(z / h, 1.0f-6)
        inner = max(1.2f0 * wstar^2 * (1.0f0 - 0.9f0*zeta) * zeta^(2.0f0/3.0f0) +
                    (1.8f0 - 1.4f0*zeta) * ustar^2, 0.0f0)
        sigw  = sqrt(inner) + 1.0f-2
        sigw  = max(sigw, 0.2f0)
        if z < aol
            tlw = 0.1f0 * z / max(sigw * (0.55f0 - 0.38f0*z/aol), 1.0f-6)
        elseif zeta < 0.1f0
            tlw = 0.59f0 * z / max(sigw, 1.0f-6)
        else
            tlw = 0.15f0 * h / max(sigw, 1.0f-6) * (1.0f0 - exp(-5.0f0*zeta))
        end
        tlw = clamp(tlw, 0.5f0, 300.0f0)
    # Stable: L > 0, pblh/|L| ≥ 1
    else
        zoh  = min(z / h, 1.0f0)
        sigw = max(1.3f0 * max(ustar, 1.0f-4) * max(1.0f0 - zoh, 0.05f0) + 1.0f-2, 0.2f0)
        tlw  = clamp(0.1f0 * h / sigw, 0.5f0, 300.0f0)
    end
    return sigw, tlw
end

function hanna_fill_tke!(metfields3D, metfields2D, indt::Integer,
                         Ni::Integer, Nj::Integer, Nk::Integer; sigw_ft::Float32 = 0.1f0)
    # Extract concrete arrays from the Dict{String,Any} ONCE: indexing the Dict
    # values inside the hot (i,j,k) loop dynamic-dispatches every access (~10^9
    # allocations).  With typed locals the loop is allocation-free.
    ustar = metfields2D["ustar"]::Array{Float32,3}
    pblh  = metfields2D["pblh"]::Array{Float32,3}
    ol    = metfields2D["ol"]::Array{Float32,3}
    zlev  = metfields3D["zlev"]::Array{Float32,4}
    tke   = metfields3D["tke"]::Array{Float32,4}
    tls   = metfields3D["tls"]::Array{Float32,4}
    it = Int(indt); Nii = Int(Ni); Njj = Int(Nj); Nkk = Int(Nk)
    @inbounds for j in 1:Njj, i in 1:Nii
        us  = ustar[i, j, it]
        pb  = pblh[i, j, it]
        olv = ol[i, j, it]
        ws  = _hanna_wstar(us, pb, olv)
        for k in 1:Nkk
            z = zlev[i, j, k, it]
            sigw, tlw = _hanna_sigw_tlw(z, pb, us, ws, olv, sigw_ft)
            tke[i, j, k, it] = 0.5f0 * sigw^2
            tls[i, j, k, it] = 0.5f0 * tlw * sigw
        end
    end
    return nothing
end

# Function barrier for the hydrostatic column integration.  Taking the arrays as
# CONCRETELY-TYPED arguments is essential: metfields3D/2D are Dict{String,Any}, and
# a `@threads` loop captures enclosing locals in a way that boxes them even when
# type-asserted — so the scalar element accesses dynamic-dispatch (billions of
# allocations).  As function args they are specialised and the loop is allocation-free.
# Hydrostatic column integration over a J-row chunk (typed args → no boxing).
function _metgrib_integrate_chunk!(jrange, zlev::Array{Float32,4}, puv::Array{Float32,4},
                             T3::Array{Float32,4}, qv3::Array{Float32,4},
                             t2m3::Array{Float32,3}, sp3::Array{Float32,3},
                             akm::Vector{Float64}, bkm::Vector{Float64},
                             akz::Vector{Float64}, bkz::Vector{Float64},
                             Ni::Int, Nk::Int, indt::Int)
   cst = 287.05/9.81
   @inbounds for j in jrange
      for i in 1:Ni
         spij  = Float64(sp3[i,j,indt])
         h     = 0.0
         tvold = Float64(t2m3[i,j,indt])
         pold  = spij
         zlev[i,j,1,indt] = 0.0f0
         for k in 1:Nk-1
            tv   = Float64(T3[i,j,k,indt]) * (1.0 + 0.608*Float64(qv3[i,j,k,indt]))
            pint = max(akm[k+1] + bkm[k+1]*spij, 0.1)   # floor p>0 (top half-level)
            dtv  = tv - tvold
            h   += abs(dtv) > 0.2 ? cst*log(pold/pint)*dtv/log(tv/tvold) :
                                    cst*log(pold/pint)*tv
            zlev[i,j,k+1,indt] = Float32(h)
            puv[i,j,k+1,indt]  = Float32(akz[k] + bkz[k]*spij)
            tvold = tv
            pold  = pint
         end
         puv[i,j,1,indt] = Float32(spij)   # surface interface pressure = sp
      end
   end
   return nothing
end

function _metgrib_integrate!(zlev, puv, T3, qv3, t2m3, sp3, akm, bkm, akz, bkz,
                             Ni::Int, Nj::Int, Nk::Int, indt::Int, nt::Int)
   if nt <= 1
      _metgrib_integrate_chunk!(1:Nj, zlev, puv, T3, qv3, t2m3, sp3, akm, bkm, akz, bkz, Ni, Nk, indt)
   else
      @sync for jr in Iterators.partition(1:Nj, cld(Nj, nt))
         Threads.@spawn _metgrib_integrate_chunk!(jr, zlev, puv, T3, qv3, t2m3, sp3,
                                                  akm, bkm, akz, bkz, Ni, Nk, indt)
      end
   end
   return nothing
end

# Run f!(krange) either serial (nt≤1) or split into nt level-chunks (typed closure → no box).
function _run_k!(nt::Int, lo::Int, hi::Int, f!::F) where {F}
   if nt <= 1 || hi < lo
      f!(lo:hi)
   else
      @sync for kr in Iterators.partition(lo:hi, max(1, cld(hi - lo + 1, nt)))
         Threads.@spawn f!(kr)
      end
   end
   return nothing
end
# Derived-field chunks (typed args; each writes its own k-slices).
function _rho_chunk!(kr, rho::Array{Float32,4}, puv::Array{Float32,4}, Tt::Array{Float32,4}, indt::Int)
   @inbounds @views @. rho[:,:,kr,indt] = puv[:,:,kr,indt] / Tt[:,:,kr,indt] / 287.05f0
   return nothing
end
function _theta_chunk!(kr, th::Array{Float32,4}, Tt::Array{Float32,4}, puv::Array{Float32,4}, indt::Int)
   @inbounds @views @. th[:,:,kr,indt] = Tt[:,:,kr,indt] * (100000.0f0 / max(puv[:,:,kr,indt], 1.0f0))^0.28573f0
   return nothing
end
# levdot [level-index/s] from the archived vertical velocity, plus drho/dz.
#
# UNITS (fixed 2026-08-12).  ERA5 param 77 is eta-dot in s^-1 — the GRIB header
# says `units = s**-1` and the magnitudes agree (|eta-dot| ~ 2e-6 s^-1 in the
# mid-troposphere = omega 0.1 Pa/s divided by dp/deta ~ 1e5 Pa).  It is NOT
# flex_extract's eta-dot*dp/deta in Pa/s, which is what `calc_etadot` produces
# and what the previous conversion assumed: dividing s^-1 by dp [Pa] yields
# s^-1 Pa^-1, not a level-index rate, and left the resolved vertical advection
# ~1e5 too small on this path.  The correct divisor is the eta THICKNESS of the
# layer, deta = d(a/p0 + b) between adjacent full levels, which is dimensionless
# and comes from the hybrid coefficients the reader already holds.
# Which convention is this file's vertical velocity in?  eta-dot is ~1e-6 s^-1
# in the mid-troposphere; eta-dot*dp/deta is ~1e-1 Pa/s — five orders apart, so
# a magnitude test is unambiguous.  Uses a high quantile rather than the mean so
# a quiescent domain cannot flip the decision, and looks only at interior levels
# (the surface/top slots are filled by extrapolation).  Threshold 1e-2 sits two
# orders clear of both populations.
# Plausible range for max|w| over the interior levels of a met field, used by
# _levdot_sanity! below.  A mid-latitude ERA5 field peaks near 0.3 m/s; deep
# convection in a global archive reaches a few m/s.  The bounds are deliberately
# wide — they are there to catch a factor ~1e5 convention error, not to police
# the meteorology.
const _LEVDOT_W_LO = 1.0f-3
const _LEVDOT_W_HI = 5.0f1

function _etadot_is_per_second(levdot::Array{Float32,4}, Nk::Int, indt::Int)
   klo = max(2, Nk ÷ 3); khi = min(Nk - 1, 2 * Nk ÷ 3)
   m = 0.0f0
   @inbounds for k in klo:khi
      v = @view levdot[:,:,k,indt]
      for x in v
         isfinite(x) || continue
         a = abs(x); a > m && (m = a)
      end
   end
   per_s = m < 1.0f-2
   @info "archived vertical velocity convention detected" max_abs=m units=(per_s ? "s^-1 (raw ECMWF eta-dot)" : "Pa/s (flex_extract eta-dot*dp/deta)")
   return per_s
end

"""
    _etadot_use_deta!(metconfig, levdot, Nk, indt) -> Bool

Cached convention decision for one met source.  The scan is O(Ni·Nj·Nk/3), so it
runs only on the FIRST file: a file set comes either from flex_extract or from
raw dissemination, never mixed.  An explicit `#info etadot_units` skips the scan
entirely.
"""
function _etadot_use_deta!(metconfig, levdot::Array{Float32,4}, Nk::Int, indt::Int)
   cached = get(metconfig, "etadot_per_second", nothing)
   cached === nothing || return cached::Bool
   v = _etadot_is_per_second(levdot, Nk, indt)
   metconfig["etadot_per_second"] = v
   return v
end

"""
    _levdot_sanity!(metconfig, levdot, zlev, Nk, indt)

Sanity check on the CONVERTED vertical velocity, one met source, once.

`_etadot_is_per_second` decides the convention from the magnitude of the
ARCHIVED field.  This checks the RESULT: the level-index rate times the local
layer thickness is a physical vertical velocity, so a convention that was
picked wrongly — or an archive that is in neither convention, or scaled by a
pre-processor — shows up here as a `w` orders of magnitude away from the
synoptic range, whatever the input looked like.  That is the case the
magnitude test cannot cover on its own: it reads an all-zero or vanishingly
small field as s^-1 by construction (see docs/ETADOT_UNITS.md).

The flex_extract / `calc_etadot` route is the one that motivates this.  It
delivers eta-dot*dp/deta in Pa/s and needs the /dp divisor, while the raw MARS
dissemination delivers eta-dot in s^-1 and needs /deta; getting it backwards
costs a factor of dp/deta ~ 1e5 in the resolved vertical transport, and the
run does not otherwise complain.
"""
function _levdot_sanity!(metconfig, levdot::Array{Float32,4}, zlev::Array{Float32,4},
                         Nk::Int, indt::Int)
   get(metconfig, "levdot_checked", false)::Bool && return nothing
   metconfig["levdot_checked"] = true
   klo = max(2, Nk ÷ 3); khi = min(Nk - 1, 2 * Nk ÷ 3)
   khi >= klo || return nothing
   wmax = 0.0f0; wsum = 0.0; nw = 0
   @inbounds for k in klo:khi, j in axes(levdot, 2), i in axes(levdot, 1)
      ld = levdot[i, j, k, indt]
      isfinite(ld) || continue
      dz = 0.5f0 * abs(zlev[i, j, k+1, indt] - zlev[i, j, k-1, indt])
      (isfinite(dz) && dz > 0.0f0) || continue
      w = abs(ld) * dz
      isfinite(w) || continue
      w > wmax && (wmax = w)
      wsum += w; nw += 1
   end
   nw == 0 && return nothing
   wmean = wsum / nw
   conv = get(metconfig, "etadot_per_second", nothing)
   convs = conv === nothing ? "unknown" :
           conv === true ? "s^-1 (raw ECMWF eta-dot, divided by deta)" :
                           "Pa/s (flex_extract eta-dot*dp/deta, divided by dp)"
   if _LEVDOT_W_LO <= wmax <= _LEVDOT_W_HI
      @info("resolved vertical velocity sanity check passed",
            max_w_m_s = wmax, mean_abs_w_m_s = wmean, convention = convs)
   else
      @warn("""resolved vertical velocity is OUTSIDE the synoptic range — the archived \
            vertical wind is probably not in the convention that was applied.  Expected \
            max|w| between $(_LEVDOT_W_LO) and $(_LEVDOT_W_HI) m/s over the interior levels \
            (a mid-latitude field peaks around 0.3 m/s); a wrong convention costs a factor \
            dp/deta ~ 1e5.  Declare the convention explicitly with `#info etadot_units` = \
            `s-1` (raw MARS/ECMWF dissemination) or `pa/s` (flex_extract / calc_etadot) in \
            the variables namelist, and see docs/ETADOT_UNITS.md.""",
            max_w_m_s = wmax, mean_abs_w_m_s = wmean, convention_applied = convs,
            levels_scanned = klo:khi)
   end
   return nothing
end

# `deta === nothing` selects the flex_extract convention (divide by dp).
function _levdrho_chunk!(kr, levdot::Array{Float32,4}, puv::Array{Float32,4},
                         rho::Array{Float32,4}, drhodz::Array{Float32,4}, zlev::Array{Float32,4},
                         deta::Union{Vector{Float64},Nothing}, indt::Int)
   @inbounds for k in kr
      lev = @view levdot[:,:,k,indt]
      if deta === nothing
         pk1 = @view puv[:,:,k+1,indt]; pk0 = @view puv[:,:,k,indt]
         @. lev /= (pk1 - pk0)
      else
         dek = Float32(deta[k])
         @. lev /= dek
      end
      dr  = @view drhodz[:,:,k,indt]
      rk1 = @view rho[:,:,k+1,indt]; rkm = @view rho[:,:,k-1,indt]
      zk1 = @view zlev[:,:,k+1,indt]; zkm = @view zlev[:,:,k-1,indt]
      @. dr = (rk1 - rkm) / (zk1 - zkm)
   end
   return nothing
end

# review 2026-09-04 A3 (item 1): unit convention of the archived surface fluxes.
# The GRIB header says stepType=accum / J m**-2 / N m**-2 s on EVERY ERA5 file,
# whether the values are raw MARS accumulations (J/m² over the step) or
# flex_extract de-accumulated per-second values (W/m², the archive on zion:
# sshf ∈ [−405, 223], |τ| ≤ 4 — physical only as W/m² and N/m²).  Dividing
# per-second values by 3600 makes |L| 3600× too large (neutral everywhere);
# taking accumulated values as instantaneous makes u* ~60× too large.  The
# test is model-agnostic (no model-name branch): every name map that carries
# sshf/ewss and no ustar/ol reaches the derivation (listvariables/ERA5.vars;
# WRF.vars / WRFv2.vars map HFX, AROME.vars sshf — instantaneous W/m² that the
# old unconditional ÷3600 also broke).  As for eta-dot, the magnitude separates
# the two populations by orders: hourly accumulated |sshf| reaches 1e5–1e6 J/m²
# (instantaneous W/m² stays < ~1.5e3, WRF HFX peaks ~600–800), hourly
# accumulated |τ| reaches 1e3–1e4 N s/m² (N/m² stays < ~10).
_sfcflux_is_accumulated(max_abs::Real, threshold::Real) = max_abs > threshold
const _SFCFLUX_SSHF_THRESHOLD = 5000.0f0     # |sshf| above → J/m² accumulated
const _SFCFLUX_TAU_THRESHOLD  = 50.0f0       # |τ|    above → N s/m² accumulated
# divisor bringing an archived flux array to per-second units
_sfcflux_divisor(max_abs::Real, threshold::Real, period_s::Real) =
   _sfcflux_is_accumulated(max_abs, threshold) ? Float32(period_s) : 1.0f0

"""
    _sfcflux_convention!(metconfig, metfields2D, indt, has_stress, has_sshf) -> (div_sshf, div_tau)

Divisors that bring sshf and ewss/nsss to per-second units (1 when the values
already are, the accumulation period [s] when they are accumulated).  Decided
ONCE per met source from the first file that carries the fields (a file set
is either flex_extract output or raw dissemination, never mixed), cached in
metconfig and logged with the numbers.  `#info sfcflux_units = per-s | accum`
skips the magnitude test.  The period comes from the GRIB header
(endStep − startStep of the sshf message, recorded by the reader); when the
header gives 0 the ERA5 hourly step (3600 s) is assumed with a warning.
"""
function _sfcflux_convention!(metconfig, metfields2D, indt::Int, has_stress::Bool, has_sshf::Bool)
   cached = get(metconfig, "sfcflux_div", nothing)
   cached === nothing || return cached::Tuple{Float32,Float32}
   _mx(a) = maximum(x -> isfinite(x) ? abs(x) : 0.0f0, a; init = 0.0f0)
   m_sshf = has_sshf   ? _mx(@view (metfields2D["sshf"]::Array{Float32,3})[:,:,indt]) : 0.0f0
   m_tau  = has_stress ? max(_mx(@view (metfields2D["ewss"]::Array{Float32,3})[:,:,indt]),
                             _mx(@view (metfields2D["nsss"]::Array{Float32,3})[:,:,indt])) : 0.0f0
   # a slot without flux data yet (analysis-time file without step-0 flux
   # messages) cannot decide: return neutral divisors WITHOUT caching
   ((has_sshf && m_sshf == 0) || (has_stress && m_tau == 0)) && return (1.0f0, 1.0f0)
   forced = get(metconfig, "sfcflux_accum", nothing)          # #info sfcflux_units
   acc_sshf = forced === nothing ? _sfcflux_is_accumulated(m_sshf, _SFCFLUX_SSHF_THRESHOLD) : forced::Bool
   acc_tau  = forced === nothing ? _sfcflux_is_accumulated(m_tau,  _SFCFLUX_TAU_THRESHOLD)  : forced::Bool
   period = 3600.0f0
   if acc_sshf || acc_tau
      step = get(metconfig, "sfcflux_step_s", nothing)
      if step isa Real && step > 0
         period = Float32(step)
      else
         @warn "surface fluxes: accumulated values but the GRIB header gives no accumulation period (endStep == startStep) — assuming 3600 s"
      end
   end
   div = (acc_sshf ? period : 1.0f0, acc_tau ? period : 1.0f0)
   has_sshf && @info "surface fluxes: sshf treated as $(acc_sshf ? "J/m² accumulated over $(round(Int, period)) s" : "W/m² (instantaneous, e.g. flex_extract ERA5 / WRF HFX / AROME — no division)") (max |sshf| = $(m_sshf), decided by $(forced === nothing ? "magnitude, model-agnostic" : "#info sfcflux_units"))"
   has_stress && @info "surface fluxes: ewss/nsss treated as $(acc_tau ? "N s/m² accumulated over $(round(Int, period)) s" : "N/m² (instantaneous — no division)") (max |τ| = $(m_tau), decided by $(forced === nothing ? "magnitude, model-agnostic" : "#info sfcflux_units"))"
   metconfig["sfcflux_div"] = div
   return div
end

function metgribfieldcalculation!(globalx::Bool,metconfig,metfields3D,metfields2D,indt::Int)
   # Performance: the old version allocated a full ~570 MB copy for every
   # `arr[:,:,:,indt]` slice plus per-level 2-D temporaries — ~17-30 GB of garbage
   # per file.  This version is allocation-free (@views + in-place + a threaded
   # scalar column sweep).  CRITICAL: metfields3D/2D are Dict{String,Any}, so the
   # array bindings are type-ASSERTED to concrete Array{Float32,N}; without that
   # the per-element loop would dynamic-dispatch every access (100-1000× slower).
   # Arithmetic is bitwise the same as before.  Concrete Vector{Float64} (the Dict
   # values are Any) so akm/bkm/akz/bkz indexing in the hot loop is type-stable.
   akm = Float64.(reverse(metconfig["ak"]))
   bkm = Float64.(reverse(metconfig["bk"]))
   akz = 0.5 .* (akm[1:end-1] .+ akm[2:end])
   bkz = 0.5 .* (bkm[1:end-1] .+ bkm[2:end])
   Ni = Int(metconfig["Ni"])
   if globalx
     Ni += 1
   end
   Nj = Int(metconfig["Nj"])
   Nk = Int(metconfig["Nk"])
   cst = 287.05/9.81
   # Worker count for the derived ("missing") 3D fields.  MET_NT=1 → serial.
   nt = clamp(parse(Int, get(ENV, "MET_NT", string(max(1, Threads.nthreads())))),
              1, max(1, Threads.nthreads()))

   # Surface pressure: ERA5 model-level dissemination provides lnsp (a 3D
   # message on level 1) instead of sp.  Derive sp = exp(lnsp) when sp was
   # not loaded directly.
   sp3  = metfields2D["sp"]::Array{Float32,3}
   if haskey(metfields2D,"lnsp") && maximum(@view sp3[:,:,indt]) <= 0.0f0
      @views sp3[:,:,indt] .= exp.((metfields2D["lnsp"]::Array{Float32,3})[:,:,indt])
   end
   t2m3 = metfields2D["2t"]::Array{Float32,3}
   t2m  = @view t2m3[:,:,indt]
   sp   = @view sp3[:,:,indt]
   T3   = metfields3D["T"]::Array{Float32,4}
   qv3  = metfields3D["qv"]::Array{Float32,4}
   zlev = metfields3D["zlev"]::Array{Float32,4}
   puv  = metfields3D["puv"]::Array{Float32,4}

   # ── Hydrostatic integration of geopotential height (zlev) and interface
   #    pressure (puv), column by column (threaded, allocation-free via the
   #    typed barrier).  Reads ORIGINAL (pre-shift) T,q; a surface level is
   #    inserted below afterwards (zlev=0, p=sp, T=2m).
   _metgrib_integrate!(zlev, puv, T3, qv3, t2m3, sp3, akm, bkm, akz, bkz, Ni, Nj, Nk, indt, nt)

   # With GRIB fields we add a surface layer to U,V,T,Q (and the CLD fields):
   # shift every level up by one (orig level k → k+1), backward in k to avoid
   # overwriting the source; each step is a contiguous in-place 2-D copy.
   for key in ("Ux","Vy","T","qv","umf","dmf","udr","ddr","tpf","cc")
      haskey(metfields3D,key) || continue
      A = metfields3D[key]::Array{Float32,4}
      for k in Nk:-1:2
         @views A[:,:,k,indt] .= A[:,:,k-1,indt]
      end
   end
   # CLD surface level: tpf keeps the lowest value (below-cloud scavenging acts
   # to the ground); mass fluxes / detrainment / cloud fraction vanish at sfc.
   for key in ("umf","dmf","udr","ddr","cc")
      haskey(metfields3D,key) && fill!(view(metfields3D[key],:,:,1,indt), 0.0f0)
   end
   haskey(metfields3D,"tpf") && (@views metfields3D["tpf"][:,:,1,indt] .= metfields3D["tpf"][:,:,2,indt])

   @views metfields3D["T"][:,:,1,indt]  .= t2m
   # Surface level (indl=1): duplicate the first real model-level wind (indl=2).
   @views metfields3D["Ux"][:,:,1,indt] .= metfields3D["Ux"][:,:,2,indt]
   @views metfields3D["Vy"][:,:,1,indt] .= metfields3D["Vy"][:,:,2,indt]
   @views metfields3D["qv"][:,:,1,indt] .= metfields3D["qv"][:,:,2,indt]

   # rho = puv / (T·Rd)   (threaded over level chunks, allocation-free)
   let rho = metfields3D["rho"]::Array{Float32,4}, Tt = metfields3D["T"]::Array{Float32,4}
      _run_k!(nt, 1, Nk, kr -> _rho_chunk!(kr, rho, puv, Tt, indt))
   end
   levdot = metfields3D["levdot"]::Array{Float32,4}
   rho    = metfields3D["rho"]::Array{Float32,4}
   drhodz = metfields3D["drhodz"]::Array{Float32,4}
   # Divisor for the archived vertical velocity — see the etadot_units note in
   # gridcheck_onefile!.  eta-dot [s^-1] divides by the dimensionless eta
   # thickness; flex_extract's eta-dot*dp/deta [Pa/s] divides by dp, which is
   # what the k-dependent `nothing` branch of _levdrho_chunk! does.
   deta = Vector{Float64}(undef, Nk)
   let p0 = 101325.0, nz = length(akz)
      for k in 1:Nk
         k0 = clamp(k, 1, nz); k1 = clamp(k + 1, 1, nz)
         d = (akz[k1] - akz[k0]) / p0 + (bkz[k1] - bkz[k0])
         deta[k] = abs(d) < 1e-12 ? -1e-12 : d
      end
   end
   use_deta = _etadot_use_deta!(metconfig, levdot, Nk, indt)
   _run_k!(nt, 2, Nk-1, kr -> _levdrho_chunk!(kr, levdot, puv, rho, drhodz, zlev,
                                              use_deta ? deta : nothing, indt))
   _levdot_sanity!(metconfig, levdot, zlev, Nk, indt)
   @views drhodz[:,:,1,indt]  .= drhodz[:,:,2,indt]
   @views drhodz[:,:,Nk,indt] .= drhodz[:,:,Nk-1,indt]

   # Surface levdot (indl=1): half of levdot[2] instead of 0, same as AROME path.
   @views levdot[:,:,1,indt] .= levdot[:,:,2,indt] .* 0.5f0
   @views zlev[:,:,Nk,indt]  .= zlev[:,:,Nk,indt] .+ 20000.0f0
   fill!(view(metfields3D["dsigwdz"],:, :, :, indt), 0.0f0)
   # Potential temperature θ = T·(p0/p)^(R/cp), p0=100000 Pa, R/cp≈0.2857 (threaded;
   # the per-cell pow makes this the heaviest broadcast).
   let th = metfields3D["theta"]::Array{Float32,4}, Tt = metfields3D["T"]::Array{Float32,4}
      _run_k!(nt, 1, Nk, kr -> _theta_chunk!(kr, th, Tt, puv, indt))
   end

   # ── Derive surface turbulence parameters from available GRIB fields ───────
   # ERA5 GRIB typically contains ewss/nsss (surface stress, N/m²) and sshf
   # (surface sensible heat flux, accumulated J/m²).  blh (boundary layer
   # height) and zust (friction velocity) may or may not be present.
   #
   # Priority: use directly-loaded fields when available; fall back to
   # deriving from stress + heat flux.
   #
   has_ustar = haskey(metfields2D, "ustar") && maximum(abs.(metfields2D["ustar"][:,:,indt])) > 0
   has_pblh  = haskey(metfields2D, "pblh")  && maximum(metfields2D["pblh"][:,:,indt])  > 0
   has_ol    = haskey(metfields2D, "ol")     && maximum(abs.(metfields2D["ol"][:,:,indt])) > 0

   has_stress = haskey(metfields2D, "ewss") && haskey(metfields2D, "nsss")
   has_sshf   = haskey(metfields2D, "sshf") && haskey(metfields2D, "2t")

   # review 2026-09-04 A3 (item 1): divisors that bring the archived fluxes to
   # per-second units (1 when the values already are W/m² / N/m², the
   # accumulation period when they are J/m² / N s/m²) — decided once per run
   # from the field magnitude, NOT from the GRIB header (see _sfcflux_convention!).
   div_sshf, div_tau = _sfcflux_convention!(metconfig, metfields2D, indt, has_stress, has_sshf)

   if !has_ustar && has_stress
       # ustar = sqrt(|τ|/ρ_sfc),  |τ| = sqrt(ewss²+nsss²) [N/m²=Pa]
       rho_sfc = metfields2D["sp"][:,:,indt] ./ metfields2D["2t"][:,:,indt] ./ 287.05f0
       tau_mag = sqrt.(metfields2D["ewss"][:,:,indt].^2 .+ metfields2D["nsss"][:,:,indt].^2) ./ div_tau
       metfields2D["ustar"][:,:,indt] .= sqrt.(max.(tau_mag ./ max.(rho_sfc, 0.001f0), 0.0f0))
       has_ustar = true
   end

   if !has_ol && has_sshf && has_ustar
       # Obukhov length L = -ρ·cp·T·u*³ / (κ·g·H)
       # sshf → W/m² through div_sshf (1 for per-second values, period for J/m²)
       # ERA5 sign: sshf < 0 for upward (unstable); L < 0 unstable
       H_flux  = metfields2D["sshf"][:,:,indt] ./ div_sshf   # W/m²; negative=upward heat
       rho_sfc = metfields2D["sp"][:,:,indt] ./ metfields2D["2t"][:,:,indt] ./ 287.05f0
       ustar3  = max.(metfields2D["ustar"][:,:,indt], 0.001f0).^3
       cp = 1004.0f0; κ = 0.4f0; g = 9.81f0
       ol_arr = -rho_sfc .* cp .* metfields2D["2t"][:,:,indt] .* ustar3 ./
                (κ .* g .* (H_flux .- sign.(H_flux) .* 1f-6))
       metfields2D["ol"][:,:,indt] .= clamp.(ol_arr, -5000.0f0, 5000.0f0)
       has_ol = true
   end

   if !has_pblh && has_ustar
       # BL height estimate: stable → Zilitinkevich (1972) h=0.4·sqrt(u*·L/f)
       # unstable/neutral → h ≈ 0.3·u*/f (Deardorff 1972), clamped to [100, 3000] m
       f_cor = 1.0f-4   # s⁻¹ (mid-latitude Coriolis approximation)
       for j in 1:Nj, i in 1:Ni
           us = Float32(metfields2D["ustar"][i,j,indt])
           ol = has_ol ? Float32(metfields2D["ol"][i,j,indt]) : 1.0f6
           if ol < 0.0f0
               # Convective: Deardorff (1972) h ~ 0.3·u*/f for neutral; use larger bound
               metfields2D["pblh"][i,j,indt] = clamp(0.3f0 * us / f_cor, 100.0f0, 3000.0f0)
           else
               # Stable: Zilitinkevich (1972) h = 0.4·sqrt(u*·L/f)
               h_stab = 0.4f0 * sqrt(us * max(ol, 1.0f0) / f_cor)
               metfields2D["pblh"][i,j,indt] = clamp(h_stab, 50.0f0, 2000.0f0)
           end
       end
       has_pblh = true
   end

   # ── Hanna pre-compute (option a) ──────────────────────────────────────────
   # Fill tke and tls with Hanna (1982) values so the Thomson/Langevin path
   # can run normally.  Falls back to zero (sigw floored at 0.2 m/s in
   # tke_partition) when no surface parameters are available.
   if has_ustar && has_pblh
       hanna_fill_tke!(metfields3D, metfields2D, indt, Ni, Nj, Nk)
   else
       fill!(view(metfields3D["tke"], :, :, :, indt), 0.0f0)
       fill!(view(metfields3D["tls"], :, :, :, indt), 0.0f0)
   end

   return nothing
end #metfieldcalculation
function metaromefieldcalculation!(globalx::Bool,metconfig,metfields3D,metfields2D,indt)
   # Broadcast-in-place rewrite (2026-07-19): same vectorised structure as the
   # original (plane-major, fused kernels — the scalar-loop rewrite hit an
   # unexplained per-element boxing and was reverted), but every temporary is a
   # cached 2-D buffer and every assignment is an in-place .= — this removes the
   # ~1 GiB/file of allocation churn (slices, boolean masks, allocated shifts).
   if !haskey(metconfig, "_akm_cache")
      akz0 = collect(Float64, reverse(Array(metconfig["ak"])))
      bkz0 = collect(Float64, reverse(Array(metconfig["bk"])))
      Nk0  = Int(metconfig["Nk"])
      akm0 = zeros(Nk0); bkm0 = zeros(Nk0)
      akm0[1] = 0.0f0; bkm0[1] = 1.0f0
      akm0[2:end-1] = 0.5 .* (akz0[1:end-1] .+ akz0[2:end])
      bkm0[2:end-1] = 0.5 .* (bkz0[1:end-1] .+ bkz0[2:end])
      akm0[end] = 1.0f0; bkm0[end] = 0.0f0
      metconfig["_akm_cache"] = (Float32.(akz0), Float32.(bkz0), Float32.(akm0), Float32.(bkm0))
   end
   akz, bkz, akm, bkm = metconfig["_akm_cache"]::NTuple{4,Vector{Float32}}

   Ni = globalx ? Int(metconfig["Ni"]) + 1 : Int(metconfig["Ni"])
   Nj = Int(metconfig["Nj"]); Nk = Int(metconfig["Nk"])

   if !haskey(metconfig, "_metarome_bufs")
      metconfig["_metarome_bufs"] = Dict{String,Matrix{Float32}}(
         k => Matrix{Float32}(undef, Ni, Nj) for k in ("sp","tv","tvo","po","pint","h2","dz2","r1","r2"))
   end
   B = metconfig["_metarome_bufs"]::Dict{String,Matrix{Float32}}
   sp = B["sp"]; tv = B["tv"]; tvo = B["tvo"]; po = B["po"]
   pint = B["pint"]; h2 = B["h2"]; dz2 = B["dz2"]

   t2m  = @view (metfields2D["2t"]::Array{Float32,3})[:,:,indt]
   lspv = @view (metfields2D["lsp"]::Array{Float32,3})[:,:,indt]
   zlev = @view (metfields3D["zlev"]::Array{Float32,4})[:,:,:,indt]
   puv  = @view (metfields3D["puv"]::Array{Float32,4})[:,:,:,indt]
   T3   = @view (metfields3D["T"]::Array{Float32,4})[:,:,:,indt]
   Q3   = @view (metfields3D["qv"]::Array{Float32,4})[:,:,:,indt]
   cst  = Float32(287.05/9.81)

   @. sp  = exp(lspv)
   @. h2  = 0.0f0
   @. tvo = t2m
   @. po  = sp
   @views zlev[:,:,1] .= 0.0f0
   # Hypsometric integration, plane by plane (reads PRE-shift T/q).
   # ifelse evaluates both branches elementwise: when tv≈tvold the (tv-tvold)/log(tv/tvold)
   # branch may produce NaN but is not selected — same result as the original mask form.
   for k in 1:Nk-1
      tk = @view T3[:,:,k]; qk = @view Q3[:,:,k]
      ak1 = akm[k+1]; bk1 = bkm[k+1]
      @. tv   = tk*(1.0f0 + 0.608f0*qk)
      @. pint = max(ak1 + bk1*sp, 0.1f0)
      @. h2  += ifelse(abs(tv - tvo) > 0.2f0,
                       cst*log(po/pint)*(tv - tvo)/log(tv/tvo),
                       cst*log(po/pint)*tv)
      @views zlev[:,:,k+1] .= h2
      akzk = akz[k]; bkzk = bkz[k]
      @views @. puv[:,:,k+1] = akzk + bkzk*sp
      @. tvo = tv
      @. po  = pint
   end

   # Shift u,v,T,q up one level — in-place plane copies, descending k
   for key in ("Ux","Vy","T","qv")
      A = @view (metfields3D[key]::Array{Float32,4})[:,:,:,indt]
      for k in Nk:-1:2
         @views A[:,:,k] .= A[:,:,k-1]
      end
   end
   @views puv[:,:,1] .= sp
   @views T3[:,:,1]  .= t2m
   Ux = @view (metfields3D["Ux"]::Array{Float32,4})[:,:,:,indt]
   Vy = @view (metfields3D["Vy"]::Array{Float32,4})[:,:,:,indt]
   @views Ux[:,:,1] .= Ux[:,:,2]
   @views Vy[:,:,1] .= Vy[:,:,2]
   @views Q3[:,:,1] .= Q3[:,:,2]

   # levdot: integrate NH_dW to interface w, then convert to level rate
   lev = @view (metfields3D["levdot"]::Array{Float32,4})[:,:,:,indt]
   NHW = @view (metfields3D["NH_dW"]::Array{Float32,4})[:,:,:,indt]
   @views lev[:,:,1] .= 0.0f0
   for k in 1:Nk-1
      @views @. lev[:,:,k+1] = lev[:,:,k] - NHW[:,:,k]/9.81f0
   end
   rho = @view (metfields3D["rho"]::Array{Float32,4})[:,:,:,indt]
   @. rho = puv/T3/287.05f0
   tke  = @view (metfields3D["tke"]::Array{Float32,4})[:,:,:,indt]
   drho = @view (metfields3D["drhodz"]::Array{Float32,4})[:,:,:,indt]
   dsig = @view (metfields3D["dsigwdz"]::Array{Float32,4})[:,:,:,indt]
   for k in 2:Nk-1
      @views @. dz2 = zlev[:,:,k+1] - zlev[:,:,k-1]
      @views @. lev[:,:,k]  = lev[:,:,k+1]/dz2/2.0f0
      @views @. drho[:,:,k] = (rho[:,:,k+1] - rho[:,:,k-1])/dz2
      @views @. dsig[:,:,k] = (max(sqrt(2.0f0*tke[:,:,k+1]), 0.2f0) -
                               max(sqrt(2.0f0*tke[:,:,k-1]), 0.2f0))/dz2
   end
   @views @. lev[:,:,1] = lev[:,:,2]*0.5f0

   theta = @view (metfields3D["theta"]::Array{Float32,4})[:,:,:,indt]
   @. theta = T3*(100000.0f0/puv)^0.28573f0
   @views drho[:,:,1] .= drho[:,:,2]
   @views @. dsig[:,:,1] = (max(sqrt(2.0f0*tke[:,:,2]), 0.2f0) -
                            max(sqrt(2.0f0*tke[:,:,1]), 0.2f0)) /
                           max(zlev[:,:,2] - zlev[:,:,1], 1.0f0)
   @views drho[:,:,Nk] .= drho[:,:,Nk-1]
   @views dsig[:,:,Nk] .= dsig[:,:,Nk-1]

   # 3-level running average of dsigwdz — rolling raw planes, in place
   r1 = B["r1"]; r2 = B["r2"]
   @views r1 .= dsig[:,:,1]
   @views r2 .= dsig[:,:,2]
   for k in 2:Nk-1
      nxt = @view dsig[:,:,k+1]
      @views @. dsig[:,:,k] = (r1 + r2 + nxt)/3.0f0
      tmp = r1; r1 = r2; r2 = tmp
      @. r2 = nxt
   end

   # ── Surface parameters: ustar, ol, pblh ──────────────────────────────────
   # AROME archives none of them (listvariables/AROME.vars leaves ustar and
   # pblh blank) and this native-e path never ran the derivation, so pblh
   # stayed identically ZERO for every AROME run: the [tracers] PBL-contact
   # tag never fired (it floors at 50 m), the trajectory CSV's pblh column was
   # 0, and the kh_ft free-troposphere gate saw everything above 50 m as free
   # troposphere.  Derived here exactly as the other native-e driver (MesoNH)
   # does it: log-law ustar, bulk-Richardson ol, and the TKE-THRESHOLD pblh —
   # the archived tke is the real e, so the threshold diagnostic is the right
   # one and the Rib profile is only the fallback.  This does NOT touch the
   # 3-D tke: no Hanna fill runs here, the native e stands.  Consumers are the
   # Hanna option-b branch, the kh_ft gate, plume injection and the transport
   # tags — none of them active in a Thomson run, so a native-e transport is
   # unchanged by this and only gains the diagnostic.
   _derive_surface_params!(metfields3D, metfields2D, indt, Ni, Nj, Nk; native_tke=true)
   return nothing
end #metaromefieldcalculation

function initoutgrid(runconfig,input,metconfig)
   outputconfig = Dict()
   Niout = parse(Int32,runconfig["OUTGRID"]["NUMXGRID"][1])
   Njout = parse(Int32,runconfig["OUTGRID"]["NUMYGRID"][1])
   outputconfig["Niout"] = Niout
   outputconfig["Njout"] = Njout
   listoutheight = [parse(Float32,val) for val in split(join(runconfig["OUTGRID"]["OUTHEIGHTS"]),",") if val != ""]
   Nkout = length(listoutheight)
   # ── Output grouping: what the 4th accumulator dimension IS ────────────────
   # The grids cost 12 B/cell per column (Float32 conc + Float64 mass), i.e.
   # ~18 MB/column on a 660×160×14 grid — per-BOX columns OOM-killed a
   # 3400-box file release at 70 GB.  In EASLEE each box releases one species,
   # so the memory-safe default is one column per SPECIES; per-box columns
   # remain available (COMMAND OUTPUTGROUP / TOML [output].group = "release"),
   # and stay the default for BACKWARD runs, where the whole point is
   # attributing sensitivities to individual receptor boxes.
   nb = Int(get(runconfig, "nbreleasebox", 0))
   # species number per release box (RELEASES_CTRL SPECNUM_REL; one value = all
   # boxes) — same convention as build_release_schedule
   boxsp = fill(Int32(1), nb)
   ctrl = haskey(input, "RELEASES_CTRL") ? input["RELEASES_CTRL"] :
          haskey(runconfig, "RELEASES_CTRL") ? runconfig["RELEASES_CTRL"] : nothing
   if ctrl !== nothing && haskey(ctrl, "SPECNUM_REL")
      vals = Int[]
      for entry in ctrl["SPECNUM_REL"], tok in split(String(entry), ",")
         isempty(strip(tok)) || push!(vals, parse(Int, strip(tok)))
      end
      if length(vals) == nb
         boxsp = Int32.(vals)
      elseif length(vals) >= 1
         boxsp = fill(Int32(vals[1]), nb)
      end
   end
   grp = haskey(input["COMMAND"], "OUTPUTGROUP") ?
         lowercase(strip(String(input["COMMAND"]["OUTPUTGROUP"][1]))) : "auto"
   grp == "box" && (grp = "release")
   grp in ("auto", "") &&
      (grp = Int(get(runconfig, "direction", 1)) == -1 ? "release" : "species")
   if input["COMMAND"]["IOUTPUTFOREACHRELEASE"][1] != "1" || nb == 0
      grp = "all"                        # legacy collapse: single column
      outclass     = fill(Int32(1), max(nb, 1))
      classspecies = Int32[length(unique(boxsp)) == 1 && nb > 0 ? boxsp[1] : 0]
   elseif grp == "species"
      spl = sort(unique(boxsp))
      outclass     = Int32[findfirst(==(s), spl) for s in boxsp]
      classspecies = Int32.(spl)
   else                                  # "release": one column per box
      outclass     = collect(Int32, 1:nb)
      classspecies = boxsp
   end
   Nc = length(classspecies)
   runconfig["Nc"] = Nc
   outputconfig["outclass"]     = outclass
   outputconfig["classspecies"] = classspecies
   outputconfig["outgroup"]     = grp
   grp == "species" && nb > Nc &&
      @info "output grouped by SPECIES ([output].group)" columns=Nc release_boxes=nb grid_MB_saved=round((nb - Nc) * 12e-6 * Int(Niout) * Int(Njout) * Nkout, digits=1)
   Nage = 1
   outputconfig["concarray"] = zeros(Float32,Niout,Njout,Nkout,Nc,Nage)
   outputconfig["massarray"] = zeros(Float64, size(outputconfig["concarray"]))
   dxout = parse(Float32,runconfig["OUTGRID"]["DXOUT"][1])
   dyout = parse(Float32,runconfig["OUTGRID"]["DYOUT"][1])
   outputconfig["dlon"] = dxout
   outputconfig["dlat"] = dyout
   lon0 = parse(Float32,runconfig["OUTGRID"]["OUTLON0"][1])
   lat0 = parse(Float32,runconfig["OUTGRID"]["OUTLAT0"][1])
   outputconfig["lon0"] = lon0
   outputconfig["lat0"] = lat0

   pih = pi / 180.
   r_earth = 6.371e6
   cosfunc(y) = cos(y * pih) * r_earth
   area = zeros(Niout,Njout)
    volume = zeros(Niout,Njout,Nkout)
    typeofgrid = metconfig[1]["projection"]
    if typeofgrid == "Lambert"
          # dxout and dyout are in meters for Lambert
          for j=1:Njout
                for i=1:Niout
                      area[i,j] = dxout * dyout
                end
          end
    else
          for j=1:Njout
             ylat1 = lat0+dyout*(j-1)
             ylat2 = lat0+dyout*j
             if ylat1 < 0.0 && ylat2 > 0.0
                  hzone = dyout * r_earth * pih
             else
                  cosfactp = cosfunc(ylat2)
                  cosfactm = cosfunc(ylat1)
                  if cosfactp < cosfactm
                        hzone = sqrt(r_earth^2 - cosfactp^2) - sqrt(r_earth^2 - cosfactm^2)
                  else
                        hzone = sqrt(r_earth^2 - cosfactm^2) - sqrt(r_earth^2 - cosfactp^2)
                  end
             end
             gridarea = 2.0 * pi * r_earth * hzone * dxout / 360.0
             for i=1:Niout
                area[i,j] = gridarea
             end
          end
    end
   # ── CRS bookkeeping for the NetCDF writer ────────────────────────────────
   # Same discriminator as the area formula above: on a Lambert nest the output
   # OUTLON0/DXOUT are PROJECTION METRES.  Sample the nest's own LON/LAT onto
   # the output grid so the file carries true geographic coordinates instead of
   # metres mislabelled "degrees_east" (which put the grid at 0N 0E).
   outputconfig["proj_m"] = (typeofgrid == "Lambert")
   if outputconfig["proj_m"] && haskey(metconfig[1], "longitudes")
      nlon = Float32.(metconfig[1]["longitudes"])
      nlat = Float32.(metconfig[1]["latitudes"])
      ndx  = Float64(metconfig[1]["dx"]);  ndy = Float64(metconfig[1]["dy"])
      Nx, Ny = size(nlon)
      lon2d = zeros(Float32, Niout, Njout); lat2d = zeros(Float32, Niout, Njout)
      @inbounds for j in 1:Njout, i in 1:Niout
         # output cell centre in nest projection metres → fractional nest index
         fx = clamp((lon0 + dxout * (i - 1)) / ndx + 1.0, 1.0, Float64(Nx))
         fy = clamp((lat0 + dyout * (j - 1)) / ndy + 1.0, 1.0, Float64(Ny))
         i0 = min(floor(Int, fx), Nx - 1); j0 = min(floor(Int, fy), Ny - 1)
         tx = fx - i0; ty = fy - j0
         lon2d[i, j] = (1-tx)*(1-ty)*nlon[i0,j0] + tx*(1-ty)*nlon[i0+1,j0] +
                       (1-tx)*ty*nlon[i0,j0+1]   + tx*ty*nlon[i0+1,j0+1]
         lat2d[i, j] = (1-tx)*(1-ty)*nlat[i0,j0] + tx*(1-ty)*nlat[i0+1,j0] +
                       (1-tx)*ty*nlat[i0,j0+1]   + tx*ty*nlat[i0+1,j0+1]
      end
      outputconfig["lon2d"] = lon2d
      outputconfig["lat2d"] = lat2d
      @info "output grid CRS = native (projection metres)" x_m=(lon0, lon0+dxout*(Niout-1)) lon_deg=(round(minimum(lon2d),digits=3), round(maximum(lon2d),digits=3)) lat_deg=(round(minimum(lat2d),digits=3), round(maximum(lat2d),digits=3))
   end
   if any(area .< 0)
      @warn "Negative values found in area" count=sum(area .< 0)
   end
   @info "Area statistics" min_area=minimum(area) max_area=maximum(area) mean_area=mean(area)
# OUTHEIGHTS values ARE the layer edge heights (the bottom/0m edge is implicit or
# listed first). outinterfaces is used by find_out_layer / test_z to bin parcels;
# it must equal outheight exactly — NOT midpoints between adjacent centers.
outinterfaces = listoutheight  # Float32 vector, already includes the 0m bottom edge

# Layer thicknesses for volume computation.
#
# outheight lists N interface heights: outheight[1] (bottom edge) … outheight[N] (top edge).
# find_out_layer now assigns z ∈ [outheight[k], outheight[k+1]) to slot k+1 (k=1..N-1).
# Slot k=1 is therefore NEVER populated (it is the "bottom edge marker" at z=outheight[1]).
# We still give it the correct dz so the assertion passes; it just never contributes to output.
for k in 1:Nkout
    dz = if k == 1
        # k=1 is the degenerate surface-marker slot (never populated by find_out_layer).
        # Give it the thickness of the first real layer above it so the assert passes.
        Nkout > 1 ? listoutheight[2] - listoutheight[1] : 1.0f0
    else
        listoutheight[k] - listoutheight[k-1]
    end
    dz <= 0 && (dz = 1.0f0)   # safety sentinel (should never trigger with valid OUTHEIGHTS)
    @assert dz > 0 "Non-positive output layer thickness at k=$k (listoutheight=$listoutheight)"
    volume[:, :, k] = area .* dz
end

outputconfig["area"] = area
outputconfig["volume"] = volume
outputconfig["outheight"] = listoutheight
outputconfig["outinterfaces"] = outinterfaces
   nest = runconfig["nest"]
   for n in 1:nest
     outputconfig[n] = Dict()
     # Wrap the met→output lon shift into (−180, 180]: a 0–360-convention met
     # crop (CAPTEX ERA5: 260–300°E) with a −180..180 OUTGRID (lon0 = −100)
     # otherwise gets a 360° shift and every deposit misses the grid.
     outputconfig[n]["xoutshift"] = mod(metconfig[n]["longitudes"][1] - lon0 + 180.0, 360.0) - 180.0
     outputconfig[n]["youtshift"] = metconfig[n]["latitudes"][1] - lat0
   end
   outputconfig["outheight"] = listoutheight

# zoutlinear maps physical height [m] → fractional layer index.
# outinterfaces has Nkout elements (the edge heights); layer k spans
# (outinterfaces[k-1], outinterfaces[k]], so we map outinterfaces → 1:Nkout.
zoutlinear = linear_interpolation(outinterfaces, 1:Nkout)
outputconfig["zoutlinear"] = zoutlinear
   outputconfig["grid"] = zeros(Float32,(Niout,Njout,Nkout,Nc))
   if gpu   # NOTE: relies on a global `gpu` flag defined by the caller
     outputconfig["Cugrid"] = adapt(CuArray,outputconfig["grid"])
   end
   outputconfig["relbox"] = runconfig["Nc"]
   return outputconfig
end

# Parse a listvariables/<NAME>.vars namelist: #3D / #2D / #info sections map
# EASLEE field names to input-file variable names; #info entries are copied
# into `input` (companions, leveltype, multistep, vmix keys, ...).
function extract_variables!(input)
    easlee_info = []
    input_info = []
    easlee_vars_3D = []
    input_vars_3D = []
    easlee_vars_2D = []
    input_vars_2D = []
    file_path = input["variablesfile"]
    current_mode = ""
    
    open(file_path, "r") do file
        for line in eachline(file)
            
            # Check for mode change before removing comments
            if startswith(strip(line), "#3D")
                current_mode = "3D"
                continue
            elseif startswith(strip(line), "#2D")
                current_mode = "2D"
                continue
            elseif startswith(strip(line), "#info")
                current_mode = "info"
                continue
            end
            # Remove comments and strip whitespace
            line = strip(split(line, "#")[1])
            
            # Skip empty lines
            if isempty(line)
                continue
            end
            
            # Split by equal sign and strip whitespace
            parts = split(line, "=")
            easlee_var = strip(parts[1])
            input_var = length(parts) < 2 ? "" : strip(parts[2])
            
            # Add to appropriate lists based on current mode
            if current_mode == "3D"
                push!(easlee_vars_3D, easlee_var)
                push!(input_vars_3D, input_var)
            elseif current_mode == "2D"
                push!(easlee_vars_2D, easlee_var)
                push!(input_vars_2D, input_var)
            elseif current_mode == "info"
                push!(easlee_info,easlee_var)
                push!(input_info, input_var)
            end
        end
    end
    convertnameeaslee = Dict()
    convertnameinput = Dict()
    for (name1,name2) in zip(input_vars_3D,easlee_vars_3D)
      convertnameeaslee[name2] = name1
    end
    for (name1,name2) in zip(input_vars_2D,easlee_vars_2D)
      convertnameeaslee[name2] = name1
    end
    for (name1,name2) in zip(input_vars_3D,easlee_vars_3D)
      convertnameinput[name1] = name2
    end
    for (name1,name2) in zip(input_vars_2D,easlee_vars_2D)
      convertnameinput[name1] = name2
    end
    easlee_vars_3D = filter(x -> x != "", easlee_vars_3D)
    easlee_vars_2D = filter(x -> x != "", easlee_vars_2D)
    input_vars_3D = filter(x -> x != "", input_vars_3D)
    input_vars_2D = filter(x -> x != "", input_vars_2D)
    for (name1,name2) in zip(easlee_info,input_info)
      input[name1] = name2
    end
    return convertnameeaslee,convertnameinput, easlee_vars_3D, input_vars_3D, easlee_vars_2D, input_vars_2D
end


# Per-model derived-field calculations (MERRA2, GFS, AROME-HDF, MesoNH)
include(joinpath(@__DIR__, "metcalc_models.jl"))

# Fix C: vertical-mixing closure selection from the variables namelist
include(joinpath(@__DIR__, "vmix_namelist.jl"))

# Asynchronous met prefetching (MetPrefetcher, readmetfiles_prefetch!)
include(joinpath(@__DIR__, "prefetch.jl"))
