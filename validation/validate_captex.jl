#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
# validate_captex.jl — CAPTEX-83 validation of EASLEE against DATEM observations.
#
# CAPTEX (Cross-Appalachian Tracer Experiment, Sep–Oct 1983): PMCH tracer
# releases from Dayton OH and Sudbury ON; surface air concentrations measured at
# ~84 receptor stations.  We run EASLEE forward for a release, sample the modelled
# surface concentration at each station over each sampling window, and compute the
# ATMES-II / DATEM statistics (FB, NMSE, R, FAC2, FMS, KSP, RANK).
#
#   OBS   directory (emit-t1.txt, meas-t1.txt, stations.txt):  ENV CAPTEX_OBS
#   MET   root (per-day ERA5 dirs IFS_EA_*/IFS_SFC_*):         ENV CAPTEX_MET=/data/CAPTEX
#
#   MODE=selftest                       parse the obs files + self-test the stats (no model)
#   MODE=run DIRECTION=forward  RELEASE=all|N|"1,3"  [RUNDAYS=4 NP=200000]
#                                       run each release from its source, sample all
#                                       stations, pool ATMES-II stats (+ per-release table)
#   MODE=run DIRECTION=backward [MAXBACK=8]
#                                       adjoint check: release backward from the strongest
#                                       receptors, convolve the SRR with the emissions
#
# NOTE — met/obs YEAR: CAPTEX is 1983; the ERA5 under /data/CAPTEX is 1984 (same
# month/day).  The script auto-detects the met year and shifts the run + obs to it
# so the pipeline runs, BUT 1984 meteorology vs 1983 tracer transport is NOT a
# valid validation — supply 1983 ERA5 for real numbers.  A loud warning is printed.
# ─────────────────────────────────────────────────────────────────────────────

using Dates, Printf, Statistics

const OBS_DIR = get(ENV, "CAPTEX_OBS",
    get(ENV, "EASLEE_CAPTEX_DIR", joinpath(@__DIR__, "..", "data", "CAPTEX")))
const MET_ROOT = get(ENV, "CAPTEX_MET", "/data/CAPTEX")
const MODE     = get(ENV, "MODE", "selftest")

# DATEM parsers + ATMES-II statistics — shared with check_etex.jl.
include(joinpath(@__DIR__, "datem.jl"))

# ── self-test: parsers on the real files + stats on synthetic data ────────────
function selftest()
    em = read_emit(joinpath(OBS_DIR, "emit-t1.txt"))
    st = read_stations(joinpath(OBS_DIR, "stations.txt"))
    ob = read_meas(joinpath(OBS_DIR, "meas-t1.txt"))
    @printf("emit:     %d releases  (first %s @ %.2f,%.2f  %.0f g)\n",
            length(em), em[1].start, em[1].lat, em[1].lon, em[1].mass_g)
    @printf("stations: %d\n", length(st))
    nz = count(o -> o.conc_pg > 0, ob)
    @printf("meas:     %d records (%d nonzero, max %.1f pg/m³, %s … %s)\n",
            length(ob), nz, maximum(o.conc_pg for o in ob), ob[1].start, ob[end].stop)
    # stats self-test: perfect, biased, and uncorrelated cases
    o = Float64[0, 1, 2, 5, 10, 20, 50, 100]
    @printf("\nstats self-test (rank↑ better, 0–4):\n")
    for (name, m) in (("perfect", copy(o)), ("×2 high", 2 .* o), ("shuffled", reverse(o)),
                      ("half-zero", [0,0,0,0,10,20,50,100.0]))
        s = datem_stats(o, m)
        @printf("  %-10s  R=%.2f FB=%+.2f NMSE=%.2f FAC2=%.2f FMS=%.0f KSP=%.0f RANK=%.2f\n",
                name, s.r, s.fb, s.nmse, s.fac2, s.fms, s.ksp, s.rank)
    end
end

if MODE == "selftest"
    selftest()
else
    include(joinpath(@__DIR__, "validate_captex_run.jl"))   # heavy model path (see file)
end
