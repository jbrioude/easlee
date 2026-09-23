
# ── flexpart_import — LEGACY-MIGRATION utility: old FLEXPART setup → RunSpec ──
# import_flexpart converts an existing FLEXPART run (pathnames + COMMAND /
# RELEASES / OUTGRID option files) into a RunSpec, ONE WAY, so users migrating
# from FLEXPART can write it back out as TOML (write_runspec) and never touch
# the legacy files again.  EASLEE runs themselves never read these files.
# A self-contained namelist parser (the FLEXPART files are simple key=value
# blocks) keeps this loadable without the Main-scope legacy reader.

module FlexpartImport

using Dates
using ..RunSpecs

export import_flexpart

# Parse &BLOCK ... / namelists into Dict(block => Dict(key => [values...])).
function _read_namelist(path::AbstractString)
    nl = Dict{String,Dict{String,Vector{String}}}()
    cur = ""; lastkey = ""
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, "&")
            cur = uppercase(strip(line[2:end])); get!(nl, cur, Dict{String,Vector{String}}())
        elseif line == "/"
            cur = ""
        elseif occursin('=', line) && occursin(r"^[A-Za-z]", line)
            k, v = split(line, '=', limit=2)
            k = uppercase(strip(k)); v = strip(rstrip(strip(v), ','))
            v = strip(v, ['"', '\'', ' '])
            get!(nl, cur, Dict{String,Vector{String}}())[k] = [String(v)]
            lastkey = k
        elseif !isempty(cur) && !isempty(lastkey) && occursin(r"^[0-9]", line)
            push!(nl[cur][lastkey], String(strip(rstrip(line, ','))))
        end
    end
    return nl
end

_first(nl, blk, key, def="") =
    (haskey(nl, blk) && haskey(nl[blk], key) && !isempty(nl[blk][key])) ? nl[blk][key][1] : def

# pathnames: 4 lines (options dir, output dir, met dir, AVAILABLE) [× nests]
function _read_pathnames(path::AbstractString)
    lines = strip.(readlines(path))
    lines = filter(!isempty, lines)
    return (options=lines[1], output=lines[2], met=lines[3], avail=lines[4])
end

# AVAILABLE line 1 = variables-namelist path → basename without extension.
function _variables_name(avail::AbstractString)
    isfile(avail) || return ""
    first_line = strip(first(Iterators.filter(!isempty, eachline(avail))))
    vf = String(split(first_line)[1])
    return splitext(basename(vf))[1]
end

"""
    import_flexpart(pathnames_file) → RunSpec

Read a FLEXPART-style run (pathnames + COMMAND/RELEASES/OUTGRID) into a RunSpec.
"""
function import_flexpart(pathnames_file::AbstractString)::RunSpec
    pn  = _read_pathnames(pathnames_file)
    opt = endswith(pn.options, "/") ? pn.options : pn.options * "/"
    base = isabspath(opt) ? opt : normpath(joinpath(dirname(pathnames_file), opt))
    cmd = _read_namelist(joinpath(base, "COMMAND"))
    rel = _read_namelist(joinpath(base, "RELEASES"))
    grd = _read_namelist(joinpath(base, "OUTGRID"))

    s = RunSpec()
    ld = tryparse(Int, _first(cmd, "COMMAND", "LDIRECT", "1"))
    s.direction = (ld === nothing || ld >= 0) ? "forward" : "backward"
    s.start = DateTime(_first(cmd,"COMMAND","IBDATE","20000101") * _first(cmd,"COMMAND","IBTIME","000000"), dateformat"yyyymmddHHMMSS")
    s.stop  = DateTime(_first(cmd,"COMMAND","IEDATE","20000101") * _first(cmd,"COMMAND","IETIME","010000"), dateformat"yyyymmddHHMMSS")
    s.timestep_s = parse(Float64, _first(cmd,"COMMAND","LSYNCTIME","900"))
    s.physics.ctl = parse(Float64, _first(cmd,"COMMAND","CTL","5"))
    s.physics.ifine = round(Int, tryparse(Float64, _first(cmd,"COMMAND","IFINE","4")))
    iout = tryparse(Int, _first(cmd,"COMMAND","IOUT","1"))
    s.output.quantity = iout in (2,10) ? "mixing_ratio" : iout in (3,11) ? "both" : "concentration"
    s.output.per_release = _first(cmd,"COMMAND","IOUTPUTFOREACHRELEASE","0") == "1"
    s.physics.convection.scheme = _first(cmd,"COMMAND","LCONVECTION","0") == "1" ? "massflux" : "off"
    s.output.interval_s = parse(Float64, _first(cmd,"COMMAND","LOUTSTEP","3600"))
    s.output.average_s  = parse(Float64, _first(cmd,"COMMAND","LOUTAVER","3600"))
    s.output.sample_s   = parse(Float64, _first(cmd,"COMMAND","LOUTSAMPLE","900"))
    s.output.path = pn.output

    # meteorology (single nest from pathnames; variables from AVAILABLE).
    # Resolve met/available paths to absolute (relative to the pathnames file's
    # directory) so the resulting config is location-independent.
    pndir = dirname(abspath(pathnames_file))
    _abs(p) = isabspath(p) ? p : normpath(joinpath(pndir, p))
    push!(s.meteorology, MetSource(
        model = "ERA5",    # assumed default — edit the written TOML if the archive differs
        path = _abs(pn.met),
        variables = _variables_name(_abs(pn.avail)),
        available = _abs(pn.avail),
    ))

    # species from SPECNUM_REL (one per box, or single applied to all)
    specnums = String[]
    if haskey(rel, "RELEASES_CTRL") && haskey(rel["RELEASES_CTRL"], "SPECNUM_REL")
        for tok in rel["RELEASES_CTRL"]["SPECNUM_REL"], t in split(tok, ",")
            isempty(strip(t)) || push!(specnums, strip(t))
        end
    end

    # release boxes — RELEASE block keys are parallel vectors
    rb = get(rel, "RELEASE", Dict{String,Vector{String}}())
    nrel = haskey(rb, "LON1") ? length(rb["LON1"]) : 0
    getv(k, i, def) = (haskey(rb, k) && i <= length(rb[k])) ? rb[k][i] : def
    seen = Set{String}()
    for i in 1:nrel
        sn = i <= length(specnums) ? specnums[i] : (isempty(specnums) ? "1" : specnums[1])
        sid = "spec$(sn)"
        if !(sid in seen)
            push!(seen, sid)
            push!(s.species, Species(id=sid, kind="gas", specnum=parse(Int, sn)))
        end
        push!(s.release, ReleaseBox(
            species = sid,
            lon = [parse(Float64, getv("LON1",i,"0")), parse(Float64, getv("LON2",i,"0"))],
            lat = [parse(Float64, getv("LAT1",i,"0")), parse(Float64, getv("LAT2",i,"0"))],
            z   = [parse(Float64, getv("Z1",i,"0")),   parse(Float64, getv("Z2",i,"0"))],
            zkind = getv("ZKIND",i,"1") == "1" ? "agl" : "asl",
            time = [DateTime(getv("IDATE1",i,"20000101")*getv("ITIME1",i,"000000"), dateformat"yyyymmddHHMMSS"),
                    DateTime(getv("IDATE2",i,"20000101")*getv("ITIME2",i,"000000"), dateformat"yyyymmddHHMMSS")],
            particles = round(Int, parse(Float64, getv("PARTS",i,"10000"))),
            mass_kg = parse(Float64, getv("MASS",i,"1")),
            comment = getv("COMMENT",i,""),
        ))
    end

    # output grid
    s.output.grid = OutputGrid(
        lon0 = parse(Float64, _first(grd,"OUTGRID","OUTLON0","-180")),
        lat0 = parse(Float64, _first(grd,"OUTGRID","OUTLAT0","-90")),
        nx = parse(Int, _first(grd,"OUTGRID","NUMXGRID","360")),
        ny = parse(Int, _first(grd,"OUTGRID","NUMYGRID","180")),
        dx = parse(Float64, _first(grd,"OUTGRID","DXOUT","1")),
        dy = parse(Float64, _first(grd,"OUTGRID","DYOUT","1")),
        heights = Float64[parse(Float64, strip(h)) for h in
                  split(join(get(get(grd,"OUTGRID",Dict()),"OUTHEIGHTS",["0,100,500,1000,3000,10000"]), ","), ",")
                  if !isempty(strip(h))],
    )
    s.name = splitext(basename(rstrip(base, '/')))[1]
    return s
end

end # module FlexpartImport
