# ── Fix C (ADAPTIVE_KERNEL §8): vertical-mixing closure from the variables
#    namelist ────────────────────────────────────────────────────────────────
# The met field named "tke" means DIFFERENT things per driver:
#   • Hanna-fill drivers (ERA5 / GFS / MERRA-2 / AROME-HDF): the field-calc
#     routines overwrite it with the Hanna fill 0.5·σ_w,H² (w-variance only).
#   • Native-e drivers (MesoNH, AROME native e): the archived full TKE
#     e = ½(σ_u²+σ_v²+σ_w²).
# This file selects, AT READ TIME, the closure that maps the archived field to
# the two quantities transport needs (σ_w² and K_z).  All closure constants
# (r_w, C_K) are NAMELIST data declared in listvariables/<NAME>.vars under
# #info, next to the variable mappings — none live in transport code:
#
#   #info
#   tke_kind = e        # archived tke is the full native TKE e
#   vmix_rw  = 0.25     # r_w : CBR vertical variance share, σ_w² = 2·e·r_w
#   vmix_ck  = 0.16     # C_K : CBR SCALAR diffusivity const, K_z = C_K·ℓ·√e
#
# The tuple returned here is stored in metconfig["vmix"] by gridcheck_onefile!
# and published to the transport kernels by build_runtime via
# EASLEE.set_vertical_mixing!.
#
# This file is deliberately standalone (no package deps) so unit tests can
# include it without the GRIB/NetCDF stack.

# Drivers whose field-calc routines OVERWRITE tke with the Hanna 0.5·σ_w² fill
# (hanna_fill_tke! via metgribfieldcalculation!/_surface_turbulence!).  A tke
# mapping in their namelists is irrelevant to transport — the fill wins.
const _VMIX_HANNA_FILL_MODELS = ("ERA5", "GFS", "MERRA2", "AROME_HDF")

"""
    vmix_from_namelist(input, tke_mapped, kz_mapped) -> (native_e::Bool, r_w::Float32, c_k::Float32)

Select the vertical-mixing closure for one met source from its parsed variables
namelist (`input` — the Dict filled by extract_variables! with the #info keys),
plus two flags derived from the variable mappings:
  tke_mapped — the namelist maps a `tke` input variable (non-empty),
  kz_mapped  — the namelist maps an archived vertical exchange coefficient
               (`kz` / `exch_h` entry, non-empty).

Errors (at read time, before any transport) when a native-TKE driver does not
declare its closure constants, and for the WRF trap below.
"""
function vmix_from_namelist(input, tke_mapped::Bool, kz_mapped::Bool)
    model = uppercase(strip(string(get(input, "model", ""))))
    kind  = lowercase(strip(string(get(input, "tke_kind", ""))))

    # ── Hanna-fill drivers / no archived tke: field ≡ 0.5·σ_w,H² ─────────────
    if model in _VMIX_HANNA_FILL_MODELS || !tke_mapped
        kind in ("", "hanna", "sigw2") || error(
            "vertical mixing: model=$model is on the Hanna-fill path (tke = 0.5·σ_w²)" *
            " but the namelist declares tke_kind=$kind — remove it or fix the model entry")
        return (false, NaN32, NaN32)
    end

    # ── WRF native TKE trap ───────────────────────────────────────────────────
    # Reconstructing K_z from a WRF TKE field is PBL-scheme dependent (the
    # closure constants differ per scheme), and the MY-family schemes archive
    # q² = 2e in the "TKE" slot — silently treating it as e halves σ_w² and
    # corrupts K_z.  The only safe source is the ARCHIVED exchange coefficient
    # EXCH_H → error at read time when it is absent.
    if model == "WRF"
        if !kz_mapped
            error("WRF native TKE: vertical mixing needs the ARCHIVED exchange " *
                  "coefficient EXCH_H (map it as `kz = EXCH_H` in the namelist); " *
                  "reconstruction from the TKE field is PBL-scheme dependent " *
                  "(MY-family archives q²=2e, not e) and is not supported")
        else
            # EXCH_H is archived but the archived-K_z closure is not wired into
            # the transport kernels yet (needs a kz met array end-to-end).
            error("WRF native TKE: archived EXCH_H detected, but the " *
                  "archived-K_z closure is not wired yet (Fix C implements the " *
                  "Hanna-fill and native-e closures only)")
        end
    end

    # ── Native-e drivers (MesoNH, AROME native e) ─────────────────────────────
    kind in ("e", "tke", "native") || error(
        "vertical mixing: model=$model maps a native tke field but the namelist " *
        "does not declare its meaning — add `tke_kind = e` plus the CBR closure " *
        "constants `vmix_rw` / `vmix_ck` to the #info block (Fix C)")

    rw = tryparse(Float32, strip(string(get(input, "vmix_rw", ""))))
    ck = tryparse(Float32, strip(string(get(input, "vmix_ck", ""))))
    (rw !== nothing && rw > 0.0f0) || error(
        "vertical mixing: model=$model with tke_kind=e needs vmix_rw > 0 " *
        "(CBR vertical variance share, σ_w² = 2·e·r_w) in the namelist #info block")
    (ck !== nothing && ck > 0.0f0) || error(
        "vertical mixing: model=$model with tke_kind=e needs vmix_ck > 0 " *
        "(CBR scalar diffusivity constant, K_z = C_K·ℓ·√e) in the namelist #info block")
    return (true, rw, ck)
end
