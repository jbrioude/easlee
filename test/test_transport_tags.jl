#!/usr/bin/env julia
# test_transport_tags.jl — transport-history tags ([tracers]):
#   increment 1 — STRATOSPHERIC residence:
#     [1] a parcel held above θ=380 K tags strat every step (strat_frac→1,
#         age_strat=0); one held below never does (frac=0, age_strat grows),
#     [2] descent clock: a parcel that leaves the stratosphere has age_strat
#         reset to 0 while there, then growing once below,
#     [3] slot hygiene: an inactive slot's tags reset; reset_tags! zeroes a slot.
#   increment 2 — PBL contact / marine-continental / convection:
#     [4] a parcel with z_agl ≤ pblh tags PBL every step (pbl_frac→1, age_pbl=0),
#         over land (lsm≥0.5) → cont_frac=1, over ocean → cont_frac=0; a parcel
#         held above the PBL never tags (pbl_frac=0, age_pbl grows),
#     [5] lofting clock: PBL for 5 steps then aloft for 5 → pbl_frac≈0.5,
#         age_pbl = 5·dt; conv_count is reset by slot hygiene / reset_tags!.
#
#   julia --project=easlee test/test_transport_tags.jl
using Printf

global gpu = false
const REPO = abspath(joinpath(@__DIR__, ".."))
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
using .EASLEE: TransportTags, apply_transport_tags!, strat_fraction,
               pbl_fraction, cont_fraction, reset_tags!

npass = 0; nfail = 0
function check(name, ok; info="")
    global npass, nfail
    ok ? (npass += 1) : (nfail += 1)
    @printf("  %-58s %s %s\n", name, ok ? "PASS" : "FAIL", info)
end

# synthetic met: 10 levels, θ(k) = 300 + 20·(k−1) → θ ≥ 380 K at levels ≥ 5;
# geometric height z(k) = 500·(k−1) m AMSL (surface at k=1 → z_agl = 500·(indl−1)).
const NZ = 10
dims = (8, 8, NZ, 2)
theta = zeros(Float32, dims); zlev = zeros(Float32, dims)
for k in 1:NZ
    theta[:, :, k, :] .= Float32(300.0 + 20.0 * (k - 1))     # 300..480 K
    zlev[:, :, k, :]  .= Float32(500.0 * (k - 1))            # 0..4500 m AMSL
end
m3 = (theta = theta, zlev = zlev)
pblh = fill(Float32(1200.0), 8, 8, 2)                        # PBL top at 1200 m everywhere
m2  = (pblh = pblh,)
LC  = EASLEE.Types.SigmaPressureLevel(Int32(NZ), zeros(Float32, NZ), zeros(Float32, NZ))
GC  = (Ni = Int32(8), Nj = Int32(8))
fctx = (nests = [(metfields3D = m3, metfields2D = m2, levelconfig = LC, gridconfig = GC)],)
# land-sea mask: land at cell (4,4), ocean everywhere else.
lsm = zeros(Float32, 8, 8); lsm[4, 4] = 1.0f0

function mkstate()
    st = EASLEE.State.ParcelState(100)
    st.status[1:3] .= 0x01
    st.nest[1:3] .= Int32(1)
    st.indx[1:3] .= 4.0f0; st.indy[1:3] .= 4.0f0
    st.indl[1] = 2.0f0        # θ=320 K → troposphere; z_agl=500 m  → PBL
    st.indl[2] = 8.0f0        # θ=440 K → stratosphere; z_agl=3500 m → free trop
    st.indl[3] = 6.0f0        # θ=400 K → stratosphere (will descend)
    return st
end

dt = 600.0
println("="^76); println("[1] stratospheric detection + fraction over 10 steps"); println("="^76)
st = mkstate(); tags = TransportTags(100; strat_theta = 380.0, lsm = lsm)
for _ in 1:10
    apply_transport_tags!(st, fctx, tags, dt, [1.0])
end
check("tropospheric parcel: strat_frac = 0, never stratospheric (age=NEVER)",
      strat_fraction(tags, 1) == 0.0f0 && tags.age_strat[1] > 1f9,
      info=@sprintf("frac %.2f age %.3g s", strat_fraction(tags, 1), tags.age_strat[1]))
check("stratospheric parcel: strat_frac = 1, age = 0",
      strat_fraction(tags, 2) == 1.0f0 && tags.age_strat[2] == 0.0f0,
      info=@sprintf("frac %.2f age %.0f", strat_fraction(tags, 2), tags.age_strat[2]))

println(); println("="^76); println("[2] descent: strat parcel drops below the θ surface"); println("="^76)
st2 = mkstate(); tags2 = TransportTags(100; strat_theta = 380.0, lsm = lsm)
for _ in 1:5; apply_transport_tags!(st2, fctx, tags2, dt, [1.0]); end     # 5 steps in strat
frac_mid = strat_fraction(tags2, 3)
st2.indl[3] = 2.0f0                                                        # descend to θ=320
for _ in 1:5; apply_transport_tags!(st2, fctx, tags2, dt, [1.0]); end     # 5 steps in trop
check("parcel strat for 5 then trop for 5 → frac ≈ 0.5",
      isapprox(strat_fraction(tags2, 3), 0.5f0; atol=1e-4),
      info=@sprintf("mid %.2f final %.2f", frac_mid, strat_fraction(tags2, 3)))
check("age_strat grows after descent (5·dt since last strat)",
      isapprox(tags2.age_strat[3], Float32(5 * dt); rtol=1e-4),
      info=@sprintf("%.0f s", tags2.age_strat[3]))

println(); println("="^76); println("[3] slot hygiene"); println("="^76)
st3 = mkstate(); tags3 = TransportTags(100; strat_theta = 380.0, lsm = lsm)
apply_transport_tags!(st3, fctx, tags3, dt, [1.0])
st3.status[2] = 0x00                                    # parcel 2 becomes inactive
apply_transport_tags!(st3, fctx, tags3, dt, [1.0])
check("freed slot reset (strat_steps=0, age=NEVER)",
      tags3.strat_steps[2] == 0 && tags3.n_steps[2] == 0 && tags3.age_strat[2] > 1f9)
reset_tags!(tags3, 1)
check("reset_tags! zeroes a slot (incl. pbl/conv)",
      tags3.strat_steps[1] == 0 && tags3.n_steps[1] == 0 && tags3.age_strat[1] > 1f9 &&
      tags3.pbl_steps[1] == 0 && tags3.land_steps[1] == 0 &&
      tags3.age_pbl[1] > 1f9 && tags3.conv_count[1] == 0)

println(); println("="^76); println("[4] PBL contact + marine/continental"); println("="^76)
# parcel 1: PBL over LAND (cell 4,4); parcel 2: above the PBL; parcel 4: PBL over OCEAN.
st4 = mkstate(); tags4 = TransportTags(100; strat_theta = 380.0, lsm = lsm)
st4.status[4] = 0x01; st4.nest[4] = Int32(1)
st4.indx[4] = 6.0f0; st4.indy[4] = 6.0f0; st4.indl[4] = 2.0f0   # ocean cell, z_agl=500 → PBL
for _ in 1:10; apply_transport_tags!(st4, fctx, tags4, dt, [1.0]); end
check("PBL parcel over land: pbl_frac=1, cont_frac=1, age_pbl=0",
      pbl_fraction(tags4, 1) == 1.0f0 && cont_fraction(tags4, 1) == 1.0f0 &&
      tags4.age_pbl[1] == 0.0f0,
      info=@sprintf("pbl %.2f cont %.2f", pbl_fraction(tags4, 1), cont_fraction(tags4, 1)))
check("free-tropospheric parcel: pbl_frac=0, age_pbl grows (never in PBL)",
      pbl_fraction(tags4, 2) == 0.0f0 && tags4.age_pbl[2] > Float32(9 * dt),
      info=@sprintf("pbl %.2f age %.0f s", pbl_fraction(tags4, 2), tags4.age_pbl[2]))
check("PBL parcel over ocean: pbl_frac=1, cont_frac=0 (marine)",
      pbl_fraction(tags4, 4) == 1.0f0 && cont_fraction(tags4, 4) == 0.0f0,
      info=@sprintf("pbl %.2f cont %.2f", pbl_fraction(tags4, 4), cont_fraction(tags4, 4)))

println(); println("="^76); println("[5] lofting clock: PBL then aloft"); println("="^76)
st5 = mkstate(); tags5 = TransportTags(100; strat_theta = 380.0, lsm = lsm)
for _ in 1:5; apply_transport_tags!(st5, fctx, tags5, dt, [1.0]); end     # 5 steps in PBL (indl=2)
st5.indl[1] = 8.0f0                                                        # loft to z_agl=3500
for _ in 1:5; apply_transport_tags!(st5, fctx, tags5, dt, [1.0]); end     # 5 steps aloft
check("PBL for 5 then aloft for 5 → pbl_frac ≈ 0.5",
      isapprox(pbl_fraction(tags5, 1), 0.5f0; atol=1e-4),
      info=@sprintf("%.2f", pbl_fraction(tags5, 1)))
check("age_pbl grows after lofting (5·dt since last PBL)",
      isapprox(tags5.age_pbl[1], Float32(5 * dt); rtol=1e-4),
      info=@sprintf("%.0f s", tags5.age_pbl[1]))

println()
@printf("%d/%d checks passed\n", npass, npass + nfail)
nfail == 0 || (println("FAIL"); exit(1))
