#!/usr/bin/env julia
# test_vertical_mixing_c.jl — unit tests for Fix C (driver-consistent vertical
# mixing, ADAPTIVE_KERNEL §8).  Synthetic arrays only — no met reads, no runs.
#
#   julia --project=easlee -t 2 test/test_vertical_mixing_c.jl
#
# Proves:
#   (a) Hanna-fill path (ERA5/GFS/MERRA-2): the dissolved tke_partition
#       (τ_Lw = K_z/σ_w², σ_w²/K_z from the reader closure) reproduces the old
#       tke_partition BIT-FOR-BIT across a range of tke/yl values including the
#       SIGW_FLOOR / TLW_FLOOR edge cases, and on a full synthetic Hanna column
#       (neutral / convective / stable regimes).
#   (b) Native-e path (MesoNH/AROME): σ_w² = 2·e·r_w and K_z = C_K·ℓ·√e with
#       (r_w, C_K) parsed from the real listvariables/MESONH.vars namelist; the
#       implied K is ~an order of magnitude below the old 2√2·ℓ·√e.
#   (c) WRF native TKE without archived EXCH_H errors at read time.

using Test

const REPO = abspath(joinpath(@__DIR__, ".."))

# Standalone includes (no GRIB/NetCDF/CUDA stack needed for these units):
include(joinpath(REPO, "easlee", "src", "types", "grid.jl"))          # Thomson/Hanna/NoTurb tags
include(joinpath(REPO, "easlee", "src", "physics", "turbulence.jl"))  # closure + tke_partition
include(joinpath(REPO, "easlee", "src", "physics", "turbulence_hanna.jl")) # hanna fill formulas
include(joinpath(REPO, "easlee", "src", "io", "vmix_namelist.jl"))    # reader closure selection

# ── Reference: the OLD tke_partition, verbatim (pre-Fix C) ────────────────────
function old_tke_partition(tke_z::Float32, yl::Float32)
    sigw = max(sqrt(2.0f0 * tke_z), SIGW_FLOOR)
    tlw  = clamp(2.0f0 * yl / sigw, TLW_FLOOR, 60.0f0)
    return tlw, sigw
end

# ── Tiny #info parser for a .vars namelist (same rules as extract_variables!) ─
function vars_info(path::AbstractString)
    info = Dict{String,String}()
    mode = ""
    for line in eachline(path)
        s = strip(line)
        startswith(s, "#3D")   && (mode = "3D";   continue)
        startswith(s, "#2D")   && (mode = "2D";   continue)
        startswith(s, "#info") && (mode = "info"; continue)
        s = strip(split(line, "#")[1])
        isempty(s) && continue
        parts = split(s, "=")
        mode == "info" || continue
        info[strip(parts[1])] = length(parts) < 2 ? "" : strip(parts[2])
    end
    return info
end

@testset "Fix C — driver-consistent vertical mixing" begin

# ══ (a) Hanna-fill path: bit-for-bit invariance ═══════════════════════════════
@testset "(a) Hanna path bit-for-bit" begin
    # Default closure is the Hanna fill
    @test vertical_mixing().native_e == false

    # Grid of values including both floors and their edges:
    #   sqrt(2·0.02) ≈ 0.2 = SIGW_FLOOR;  2·yl/σ_w = 0.5 = TLW_FLOOR near
    #   yl = 0.05 when σ_w = 0.2;  large yl hits the 60 s clamp.
    tkes = Float32[0.0, 1f-8, 1f-4, 0.0199, 0.02, 0.0201, 0.1, 0.5, 2.0, 10.0, 100.0]
    yls  = Float32[0.0, 1f-3, 0.0499, 0.05, 0.0501, 0.1, 1.0, 5.0, 50.0, 500.0, 5000.0]
    n_checked = 0
    for tke in tkes, yl in yls
        tlw_new, sigw_new = tke_partition(tke, yl)                # dissolved path
        tlw_ref, sigw_ref = old_tke_partition(tke, yl)            # historic path
        @test tlw_new  === tlw_ref                                # bitwise
        @test sigw_new === sigw_ref
        n_checked += 1
    end
    @info "Hanna grid: $n_checked (tke, yl) pairs bit-for-bit identical"

    # Explicit 3-arg (closure) form agrees with the 2-arg form
    vm = vertical_mixing()
    for tke in tkes, yl in yls
        @test tke_partition(vm, tke, yl) === tke_partition(tke, yl)
    end

    # Full synthetic Hanna-fill column, all three stability regimes.
    # The reader fill writes tke = 0.5·σ_w,H² (hanna_fill_tke!): check the
    # partition recovers Hanna's σ_w to Float32 round-off and matches the old
    # partition bit-for-bit level by level.
    z_prof = Float32.(collect(0.0:50.0:3000.0))
    nz     = length(z_prof)
    theta_prof = 300.0f0 .+ 0.004f0 .* z_prof              # weakly stable θ(z)
    regimes = [                    # (name, ustar, pblh, ol)
        ("neutral",    0.4f0, 1200.0f0,  1.0f6),
        ("convective", 0.3f0, 1500.0f0, -50.0f0),
        ("stable",     0.2f0,  300.0f0,  80.0f0),
    ]
    for (name, ust, pblh, ol) in regimes
        wst = hanna_wstar(ust, pblh, ol)
        tkecol = Vector{Float32}(undef, nz)
        sigwH  = Vector{Float32}(undef, nz)
        for k in 1:nz
            sigw, _, _ = hanna_sigw_tlw(z_prof[k], pblh, ust, wst, ol, 0.1f0)
            sigwH[k]  = sigw
            tkecol[k] = 0.5f0 * sigw^2                     # the reader fill
        end
        for k in 1:nz-1
            yl = turbulent_length_scale_v2(k, k + 1, tkecol[k], z_prof[k],
                                           z_prof, theta_prof)
            tlw_new, sigw_new = tke_partition(tkecol[k], yl)
            tlw_ref, sigw_ref = old_tke_partition(tkecol[k], yl)
            @test tlw_new === tlw_ref && sigw_new === sigw_ref
            # σ_w round-trip: fill 0.5σ² then √(2·) — Float32 round-off only
            # (where above the σ_w floor; at the floor the clamp wins, as before)
            if sigwH[k] > SIGW_FLOOR
                @test isapprox(sigw_new, sigwH[k]; rtol = 4 * eps(Float32))
            else
                @test sigw_new === SIGW_FLOOR
            end
        end
        @info "Hanna column ($name): bit-for-bit vs old partition, σ_w round-trip OK"
    end
end

# ══ (b) Native-e path: CBR closure, constants from the namelist ═══════════════
@testset "(b) native-e closure (MesoNH namelist)" begin
    # Constants must come from the shipped namelist — none in transport code.
    info = vars_info(joinpath(REPO, "listvariables", "MESONH.vars"))
    @test info["tke_kind"] == "e"
    vmix_tuple = vmix_from_namelist(info, true, false)
    @test vmix_tuple[1] === true
    r_w = vmix_tuple[2]
    c_k = vmix_tuple[3]
    @test 0.0f0 < r_w < 1.0f0        # vertical share of 2e
    @test 0.14f0 <= c_k <= 0.2f0     # CBR scalar constant range (§8 table)

    set_vertical_mixing!(true; r_w = r_w, c_k = c_k)
    @test vertical_mixing().native_e

    # Synthetic native-e column: e(z) decaying, weakly stable θ(z).
    z_prof     = Float32.(collect(10.0:50.0:2000.0))
    theta_prof = 300.0f0 .+ 0.003f0 .* z_prof
    e_prof     = 2.0f0 .* exp.(-z_prof ./ 800.0f0)     # e ∈ [0.16, 2.0] m²/s²
    nz         = length(z_prof)

    ratios = Float64[]
    for k in 1:nz-1
        e  = e_prof[k]
        yl = turbulent_length_scale_v2(k, k + 1, e, z_prof[k], z_prof, theta_prof)
        tlw, sigw = tke_partition(e, yl)

        # σ_w² = 2·e·r_w (CBR vertical partition) — up to Float32 round-off
        sigw2_ref = 2.0f0 * r_w * e
        if sqrt(sigw2_ref) > SIGW_FLOOR
            @test isapprox(sigw^2, sigw2_ref; rtol = 4 * eps(Float32))
        else
            @test sigw === SIGW_FLOOR
        end

        # K_z = σ_w²·τ_Lw = C_K·ℓ·√e when the τ_Lw clamps are inactive
        K_new     = Float64(sigw)^2 * Float64(tlw)
        K_target  = Float64(c_k) * Float64(yl) * sqrt(Float64(e))
        tlw_raw   = K_target / Float64(sigw)^2
        if TLW_FLOOR < tlw_raw < 60.0
            @test isapprox(K_new, K_target; rtol = 1e-5)

            # Old implied diffusivity: σ_old = √(2e), τ_old = 2ℓ/σ_old
            #   K_old = σ_old²·τ_old = 2·ℓ·√(2e) = 2√2·ℓ·√e
            tlw_o, sigw_o = old_tke_partition(e, yl)
            if TLW_FLOOR < 2.0f0 * yl / sigw_o < 60.0f0
                K_old = Float64(sigw_o)^2 * Float64(tlw_o)
                push!(ratios, K_old / K_new)
                # analytic ratio 2√2/C_K, independent of r_w
                @test isapprox(K_old / K_new, 2 * sqrt(2.0) / c_k; rtol = 1e-4)
            end
        end
    end
    @test !isempty(ratios)
    # "~an order of magnitude below": 2√2/0.16 ≈ 17.7
    @test all(r -> 10.0 < r < 40.0, ratios)
    @info "native-e: σ_w²=2e·r_w and K=C_K·ℓ√e verified; K_old/K_new ≈ " *
          string(round(sum(ratios)/length(ratios); digits = 2)) *
          " (2√2/C_K = $(round(2*sqrt(2.0)/c_k; digits = 2)))"

    # Floor edge: e = 0 → σ_w at floor, τ_Lw at floor (K_z → 0)
    tlw0, sigw0 = tke_partition(0.0f0, 100.0f0)
    @test sigw0 === SIGW_FLOOR && tlw0 === TLW_FLOOR

    # native-e constants are namelist-mandatory
    @test_throws ErrorException set_vertical_mixing!(true)
    @test_throws ErrorException set_vertical_mixing!(true; r_w = 0.25f0)

    # reset to the Hanna default for anything running after this test
    set_vertical_mixing!(false)
    @test !vertical_mixing().native_e
end

# ══ (c) Read-time guards (vmix_from_namelist) ═════════════════════════════════
@testset "(c) reader guards" begin
    # WRF native TKE without archived EXCH_H → error (the §8 trap: K_z
    # reconstruction is PBL-scheme dependent; MY-family archives q²=2e).
    wrf = Dict("model" => "WRF", "format" => "NETCDF")
    err = try vmix_from_namelist(wrf, true, false); nothing catch e; e end
    @test err isa ErrorException
    @test occursin("EXCH_H", err.msg)

    # WRF with EXCH_H mapped: archived-K_z closure not wired → explicit error
    @test_throws ErrorException vmix_from_namelist(wrf, true, true)

    # WRF with NO tke mapped (current WRF.vars: `tke =` empty) → Hanna path, OK
    @test vmix_from_namelist(wrf, false, false) === (false, NaN32, NaN32)

    # Native driver mapping tke but missing declaration / constants → error
    mnh = Dict("model" => "MESONH")
    @test_throws ErrorException vmix_from_namelist(mnh, true, false)  # no tke_kind
    mnh["tke_kind"] = "e"
    @test_throws ErrorException vmix_from_namelist(mnh, true, false)  # no vmix_rw
    mnh["vmix_rw"] = "0.25"
    @test_throws ErrorException vmix_from_namelist(mnh, true, false)  # no vmix_ck
    mnh["vmix_ck"] = "0.16"
    @test vmix_from_namelist(mnh, true, false) === (true, 0.25f0, 0.16f0)

    # Hanna-fill drivers: closure is (false, NaN, NaN) whatever tke says
    for m in ("ERA5", "GFS", "MERRA2", "AROME_HDF")
        @test vmix_from_namelist(Dict("model" => m), true, false)[1] === false
    end
    # ... but declaring tke_kind=e on a Hanna-fill driver is a config error
    @test_throws ErrorException vmix_from_namelist(
        Dict("model" => "ERA5", "tke_kind" => "e"), true, false)

    # The shipped AROME.vars (native e) parses to a valid closure
    ainfo = vars_info(joinpath(REPO, "listvariables", "AROME.vars"))
    @test vmix_from_namelist(ainfo, true, false)[1] === true
end

end # outer testset

println("\nAll Fix C vertical-mixing tests passed.")
