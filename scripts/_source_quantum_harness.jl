# _source_quantum_harness.jl — synthetic 1° grid with cell- and time-varying surface fields:
# run the sea-spray and dust sources for `nsteps` and compare the mass carried by the spawned
# parcels with the flux integrated INDEPENDENTLY from the pure flux functions over the same
# cells and steps.  Included (after the model is loaded) by test/test_source_mass_quantum.jl
# and scripts/paper_fig_sources_data.jl.
function source_quantum_harness(; nsteps::Int = 12, dt_s::Float64 = 900.0, seed::Int = 1,
                                  q_spray::Float64 = 5.0e3, q_dust::Float64 = 5.0e7)
    O = EASLEE.Operators
    G = EASLEE.Types.Rectilinear(Int32(360), Int32(180), Int32(10), 1.0f0, 1.0f0,
            0.0f0, -90.0f0, Float32(1/111320.0), Float32(1/111320.0),
            1.0f0, 0.0f0, 1.0f0, 0.0f0, false, false, false, 1.0f0, 1.0f0, true)
    LC = EASLEE.Types.SigmaPressureLevel(Int32(10), zeros(Float32, 10), zeros(Float32, 10))
    fctx = (nests = [(gridconfig = G, levelconfig = LC)],)
    rng = MersenneTwister(seed)
    lsm = zeros(Float32, 360, 180, 2); lsm[:, 100:130, :] .= 1.0f0            # land band ≈ 10–40 °N
    u10 = 4.0f0 .+ 12.0f0 .* rand(rng, Float32, 360, 180, 2)                # 4–16 m/s per cell and time slot
    ust = 0.15f0 .+ 0.55f0 .* rand(rng, Float32, 360, 180, 2)               # 0.15–0.70 m/s
    T2 = 285.0; PS = 98000.0
    m2 = Dict{String,Any}("lsm" => lsm, "10u" => u10, "10v" => zeros(Float32, 360, 180, 2),
                          "2t" => fill(Float32(T2), 360, 180, 2), "ustar" => ust,
                          "sp" => fill(Float32(PS), 360, 180, 2))
    t1 = DateTime(2025, 5, 1); t2 = DateTime(2025, 5, 2)
    rs = EasleeConfig.ReleaseBox(species="ss", kind="seaspray", lon=[10.0, 20.0], lat=[-20.0, 0.0],
                                 time=[t1, t2], mass_per_parcel_kg=q_spray)
    rd = EasleeConfig.ReleaseBox(species="dust", kind="dust", lon=[10.0, 20.0], lat=[12.0, 30.0],
                                 time=[t1, t2], mass_per_parcel_kg=q_dust,
                                 dust_clay=0.2, dust_bare=1.0, dust_w_pct=0.5)
    ssb = O.build_seaspray(fctx, m2, 1, rs); db = O.build_dust(fctx, m2, 2, rd)
    st = EASLEE.State.ParcelState(200_000)
    for k in 1:nsteps
        indt = isodd(k) ? 1.0 : 2.0
        O.apply_seaspray!(st, fctx, m2, ssb, t1 + Minute(15k), indt, dt_s; rng=rng)
        O.apply_dust!(st, fctx, m2, db, t1 + Minute(15k), indt, dt_s; rng=rng)
    end
    Is = 0.0; Id = 0.0; ra = PS / (287.05 * T2)
    for k in 1:nsteps
        it = isodd(k) ? 1 : 2
        for c in eachindex(ssb.ii)
            i = Int(ssb.ii[c]); j = Int(ssb.jj[c])
            Is += O.ossa_mass_flux(O.ossa_reynolds(Float64(u10[i, j, it]), T2); rho_dry=2017.0) * ssb.area[c] * dt_s
        end
        for c in eachindex(db.ii)
            i = Int(db.ii[c]); j = Int(db.jj[c])
            Id += O.dust_flux_kok(Float64(ust[i, j, it]), ra; clay=0.2, bare=1.0, w_pct=0.5) * db.area[c] * dt_s
        end
    end
    act = [i for i in eachindex(st.status) if st.status[i] == 0x01]
    ps = [i for i in act if st.relbox[i] == 1]; pd = [i for i in act if st.relbox[i] == 2]
    summ(v) = sum(Float64(st.mass[i]) for i in v; init=0.0)
    return (spray = (integral=Is, spawned=summ(ps), emitted=ssb.emitted[], unspawned=sum(ssb.acc; init=0.0),
                     n=ssb.nspawn[], nparcels=length(ps), ncells=length(ssb.ii), q=ssb.quantum, maxacc=maximum(ssb.acc; init=0.0),
                     indl=(isempty(ps) ? NaN : maximum(Float64(st.indl[i]) for i in ps))),
            dust  = (integral=Id, spawned=summ(pd), emitted=db.emitted[], unspawned=sum(db.acc; init=0.0),
                     n=db.nspawn[], nparcels=length(pd), ncells=length(db.ii), q=db.quantum, maxacc=maximum(db.acc; init=0.0),
                     indl=(isempty(pd) ? NaN : maximum(Float64(st.indl[i]) for i in pd))))
end
