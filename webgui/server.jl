using HTTP, JSON3, TOML
import Contour as _CTR

const FRONTEND_DIR = joinpath(@__DIR__, "frontend")

# ── Helpers ───────────────────────────────────────────────────────────────────
json_ok(data=Dict{String,Any}()) =
    HTTP.Response(200, ["Content-Type" => "application/json"], JSON3.write(data))

json_err(msg::String; status=400) =
    HTTP.Response(status, ["Content-Type" => "application/json"],
                  JSON3.write(Dict("ok" => false, "msg" => msg)))

const MIME_MAP = Dict(
    ".html" => "text/html; charset=utf-8",
    ".js"   => "application/javascript",
    ".css"  => "text/css",
    ".png"  => "image/png",
    ".svg"  => "image/svg+xml",
    ".ico"  => "image/x-icon",
    ".json" => "application/json",
)

# ── Route handlers ────────────────────────────────────────────────────────────
function handle_index(req::HTTP.Request)
    HTTP.Response(200, ["Content-Type" => "text/html; charset=utf-8"],
                  read(joinpath(FRONTEND_DIR, "index.html")))
end

function serve_file(req::HTTP.Request, strip_prefix::String="")
    path = HTTP.URI(req.target).path
    isempty(strip_prefix) || (path = replace(path, Regex("^" * strip_prefix) => ""))
    path = replace(path, r"\.\." => "")          # block path traversal
    path = lstrip(path, '/')
    file = joinpath(FRONTEND_DIR, path)
    isfile(file) || return HTTP.Response(404, "Not found")
    ext  = lowercase(splitext(file)[2])
    mime = get(MIME_MAP, ext, "application/octet-stream")
    HTTP.Response(200, ["Content-Type" => mime], read(file))
end

# /static/foo  →  frontend/foo   (legacy; kept for compatibility)
handle_static(req::HTTP.Request)  = serve_file(req, "/static")
# /app.js, /style.css  →  frontend/app.js, frontend/style.css
handle_asset(req::HTTP.Request)   = serve_file(req)

function handle_state(::HTTP.Request)
    json_ok(Dict(
        "running"  => STATE["running"],
        "step"     => STATE["step"],
        "time"     => STATE["time"],
        "progress" => STATE["progress"],
        "nactive"  => STATE["nactive"],
        "log"      => get_log_snapshot(),
    ))
end

function handle_figure(req::HTTP.Request)
    p        = HTTP.queryparams(HTTP.URI(req.target))
    view     = Symbol(get(p, "view",      "column"))
    field    = get(p, "met_field", "Ux")
    level    = parse(Int,     get(p, "level",    "1"))
    lat_frac = parse(Float64, get(p, "lat_frac", "0.5"))
    sp       = get(p, "species", "");  species = isempty(sp) ? nothing : sp
    # Reload the requested view on demand (cheap if already cached)
    view == :met && load_met_from_live!(field, level)
    view == :column && !isempty(LIVE.last_nc) &&
        load_column_from_nc!(LIVE.last_nc; species=species)
    view == :cross && !isempty(LIVE.last_nc) &&
        load_crosssection_from_nc!(LIVE.last_nc, lat_frac; species=species)
    png = render_figure(view; met_field=field, level=level, lat_frac=lat_frac)
    HTTP.Response(200, ["Content-Type" => "image/png",
                        "Cache-Control" => "no-store"], png)
end

function handle_control(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    body   = JSON3.read(String(req.body))
    action = String(body["action"])
    if action == "run"
        STATE["running"] &&
            return json_err("already running")
        while isready(CONTROL_CH); take!(CONTROL_CH); end
        SPEC_DICT[] === nothing && (SPEC_DICT[] = default_spec_dict())
        cfg = deepcopy(SPEC_DICT[])
        Threads.@spawn run_model_spec!(cfg)
        return json_ok(Dict("ok" => true))
    elseif action == "pause"
        pause_model!()
        return json_ok(Dict("ok" => true))
    elseif action == "resume"
        resume_model!()
        return json_ok(Dict("ok" => true))
    elseif action == "stop"
        stop_model!()
        return json_ok(Dict("ok" => true))
    end
    json_err("unknown action: $action")
end

function handle_eval(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    get(ENV, "EASLEE_GUI_REPL", "0") == "1" ||
        return json_err("the Julia console is disabled; start the GUI with EASLEE_GUI_REPL=1 to enable it"; status=403)
    body   = JSON3.read(String(req.body))
    code   = String(body["code"])
    out, is_err = eval_expression(code)
    json_ok(Dict("output" => out, "error" => is_err))
end

# GET /api/browse?path=DIR → list a server-side directory for the file/folder
# picker dialog.  Local single-user tool, so browsing the filesystem is the
# point (no sandbox).  Returns dirs-first entries plus the parent path.
function handle_browse(req::HTTP.Request)
    q    = HTTP.queryparams(HTTP.URI(req.target))
    path = get(q, "path", "")
    isempty(path) && (path = homedir())
    path = abspath(expanduser(path))
    isdir(path) || (path = dirname(path))            # a file → show its folder
    isdir(path) || return json_err("not a directory: $path"; status=404)
    entries = Dict{String,Any}[]
    try
        for name in readdir(path)
            startswith(name, ".") && continue        # hide dotfiles
            push!(entries, Dict("name" => name, "dir" => isdir(joinpath(path, name))))
        end
    catch e
        return json_err("cannot read $path: $e"; status=403)
    end
    sort!(entries, by = e -> (!e["dir"], lowercase(e["name"])))   # dirs first, then files
    parent = dirname(path)
    json_ok(Dict("path" => path, "parent" => (parent == path ? "" : parent),
                 "sep" => string(Base.Filesystem.path_separator), "entries" => entries))
end

# GET /api/domain → lon/lat outline + bounding box of each meteorology nest in
# the current config, for the release-placement map.  Reuses build_run for the
# projection (same path a run uses), cached per meteorology block so it builds at
# most once per distinct met config; reuses a live runtime when one exists.
# Domain cache lives IN the met directory (.easlee_domain.json), so it survives
# restarts and the expensive grid build runs at most once per dataset.  One file
# can hold several datasets sharing the directory, keyed by model+variables.
const _DOMAIN_FILE = ".easlee_domain.json"
_domain_key(m) = "$(m.model)|$(m.variables)"

function _load_domain_cache(dir::AbstractString)
    f = joinpath(dir, _DOMAIN_FILE)
    isfile(f) || return Dict{String,Any}()
    try
        return copy(JSON3.read(read(f, String), Dict{String,Any}))
    catch
        return Dict{String,Any}()
    end
end

function _save_domain_entry(dir::AbstractString, key::AbstractString, payload)
    try
        isdir(dir) || return
        cache = _load_domain_cache(dir)
        cache[key] = payload
        open(joinpath(dir, _DOMAIN_FILE), "w") do io
            JSON3.write(io, cache)
        end
    catch
        # read-only met directory etc. — caching is best-effort, the map still works
    end
end

function _nest_outline(gc)
    Ni = Int(gc.Ni); Nj = Int(gc.Nj); np = 40
    pts = Tuple{Float64,Float64}[]   # (lon, lat)
    addpt(i, j) = (latlon = EASLEE.Operators.indyindxtolatlon(gc, Float64(j), Float64(i));
                   push!(pts, (Float64(latlon[2]), Float64(latlon[1]))))
    is = range(1, Ni; length=np); js = range(1, Nj; length=np)
    for i in is; addpt(i, 1);  end          # south edge
    for j in js; addpt(Ni, j); end          # east edge
    for i in reverse(is); addpt(i, Nj); end # north edge
    for j in reverse(js); addpt(1, j);  end # west edge
    lons = [p[1] for p in pts]; lats = [p[2] for p in pts]
    glob = hasproperty(gc, :globalx) && gc.globalx
    lonr = glob ? [-180.0, 180.0] : [minimum(lons), maximum(lons)]
    latr = glob ? [ -90.0,  90.0] : [minimum(lats), maximum(lats)]
    return Dict("ni"=>Ni, "nj"=>Nj, "global"=>glob,
                "outline"=>[[p[1], p[2]] for p in pts],
                "lon_range"=>lonr, "lat_range"=>latr)
end

# Coastline polylines for the release map: the 0.5 contour of the met LAND-SEA
# MASK — met-grid resolution, i.e. far better than the bundled low-res
# coastline.json.  Needs one met time step read (lsm is a read field); returns
# nothing when the dataset has no lsm or the read fails.
function _lsm_coast(rt, spec, n)
    try
        m2 = rt.metfields2D[n]
        haskey(m2, "lsm") || return nothing
        lsm = m2["lsm"]
        maximum(x -> isfinite(x) ? abs(x) : 0.0f0, view(lsm, :, :, 1)) > 0 ||
            read_met_step!(rt, spec.start)                # fill the arrays once
        gc = rt.ctx.nests[n].gridconfig
        Ni = min(Int(gc.Ni), size(lsm, 1)); Nj = min(Int(gc.Nj), size(lsm, 2))
        lons = Float64[EASLEE.Operators.indyindxtolatlon(gc, 1.0, Float64(i))[2] for i in 1:Ni]
        lats = Float64[EASLEE.Operators.indyindxtolatlon(gc, Float64(j), 1.0)[1] for j in 1:Nj]
        (issorted(lons) && issorted(lats)) || return nothing
        c = _CTR.contour(lons, lats, Float64.(lsm[1:Ni, 1:Nj, 1]), 0.5)
        coast = Vector{Vector{Vector{Float64}}}()
        for ln in _CTR.lines(c)
            xs, ys = _CTR.coordinates(ln)
            length(xs) < 2 && continue
            stp = max(1, length(xs) ÷ 600)                # cap points per line
            push!(coast, [[round(xs[i]; digits=3), round(ys[i]; digits=3)]
                          for i in 1:stp:length(xs)])
        end
        return isempty(coast) ? nothing : coast
    catch e
        @warn "lsm coastline failed: $(typeof(e))"
        return nothing
    end
end

function handle_domain(::HTTP.Request)
    SPEC_DICT[] === nothing && (SPEC_DICT[] = default_spec_dict())
    local spec
    try
        spec = dict_to_runspec(SPEC_DICT[])
    catch e
        return json_err("invalid configuration: $e")
    end
    isempty(spec.meteorology) && return json_err("no meteorology in configuration")
    nsrc  = length(spec.meteorology)
    nests = Vector{Any}(undef, nsrc)
    miss  = false
    for (n, m) in enumerate(spec.meteorology)             # try the per-directory cache first
        cache = _load_domain_cache(m.path)
        k = _domain_key(m)
        # entries cached before the lsm-coastline feature lack "coast" → rebuild once
        if haskey(cache, k) && haskey(cache[k], "coast")
            nests[n] = cache[k]
        else
            miss = true
        end
    end
    if miss                                               # build once, then cache to each met dir
        try
            rt = first(build_run(spec; prefetch=false))
            for n in 1:min(nsrc, length(rt.ctx.nests))
                m = spec.meteorology[n]
                payload = merge(_nest_outline(rt.ctx.nests[n].gridconfig), Dict("model"=>string(m.model)))
                coast = _lsm_coast(rt, spec, n)
                payload["coast"] = coast === nothing ? Any[] : coast
                nests[n] = payload
                _save_domain_entry(m.path, _domain_key(m), payload)
            end
        catch e
            return json_err("could not read met domain: $e"; status=500)
        end
    end
    HTTP.Response(200, ["Content-Type"=>"application/json"], JSON3.write(Dict("nests"=>collect(nests))))
end

# ── Schema-driven config endpoints ────────────────────────────────────────────
# /api/schema → the field descriptor list (drives the auto-generated form);
# /api/config → the current RunSpec as a TOML-shaped dict (GET), or set it (POST).
# the Chemistry section is not offered: the reactive-chemistry module is not part of this distribution
_schema_json_transport() = filter(f -> get(f, "group", get(f, :group, "")) != "Chemistry", schema_json())
handle_schema(::HTTP.Request) =
    HTTP.Response(200, ["Content-Type"=>"application/json"], JSON3.write(_schema_json_transport()))

function handle_get_config(::HTTP.Request)
    SPEC_DICT[] === nothing && (SPEC_DICT[] = default_spec_dict())
    HTTP.Response(200, ["Content-Type"=>"application/json"], JSON3.write(SPEC_DICT[]))
end

function handle_set_config(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    try
        cfg = JSON3.read(String(req.body), Dict{String,Any})
        validate_runspec(dict_to_runspec(cfg))    # validate before storing
        SPEC_DICT[] = cfg
        json_ok(Dict("ok"=>true))
    catch e
        json_err("invalid configuration: $e")
    end
end

# POST /api/load_toml {path} → parse+validate a .toml config from disk, store it
# as the active RunSpec, and return it as the TOML-shaped dict the form renders.
function handle_load_toml(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    body = JSON3.read(String(req.body))
    path = String(body["path"])
    isfile(path) || return json_err("file not found: $path"; status=404)
    try
        spec = read_runspec(path)                  # parse + unit-convert + validate
        d = runspec_to_dict(spec)
        SPEC_DICT[] = d
        HTTP.Response(200, ["Content-Type"=>"application/json"], JSON3.write(d))
    catch e
        json_err("TOML load error: $e"; status=500)
    end
end

# GET /api/config.toml → the current configuration serialised as a TOML document
# (normalised through RunSpec so it is exactly what `read_runspec` would accept).
function handle_config_toml(::HTTP.Request)
    SPEC_DICT[] === nothing && (SPEC_DICT[] = default_spec_dict())
    try
        spec = dict_to_runspec(SPEC_DICT[])
        io = IOBuffer()
        println(io, "# EASLEE run configuration")
        TOML.print(io, runspec_to_dict(spec); sorted=true)
        HTTP.Response(200, ["Content-Type"=>"text/plain; charset=utf-8"], take!(io))
    catch e
        HTTP.Response(500, "serialise error: $e")
    end
end

# ── Router ────────────────────────────────────────────────────────────────────

# ── security (review 2026-09-23) ─────────────────────────────────────────────
# A local single-user tool: bind to the loopback interface unless EASLEE_GUI_HOST
# says otherwise; refuse POSTs that are not JSON from a local Host (a page in the
# browser cannot send a cross-site JSON POST without a preflight, and the Host
# header must be local); the Julia console (/api/eval) is OFF unless EASLEE_GUI_REPL=1.
const _LOCAL_HOSTS = ("localhost", "127.0.0.1", "[::1]", "::1")
function _host_name(h::AbstractString)
    h = strip(h)
    startswith(h, "[") && return h[1:findfirst(==(']'), h)]      # [::1]:port
    return first(split(h, ':'))
end
function _post_guard(req::HTTP.Request)
    ct = lowercase(HTTP.header(req, "Content-Type", ""))
    startswith(ct, "application/json") ||
        return json_err("POST requests must carry Content-Type: application/json"; status=415)
    host = _host_name(HTTP.header(req, "Host", ""))
    host in _LOCAL_HOSTS ||
        return json_err("requests are accepted from localhost only (Host = $host)"; status=403)
    return nothing
end

# ── Meteorology download (download/*.jl) ──────────────────────────────────────
# POST /api/met/download  {model: GFS|MERRA2|ERA5, out, ...}  → spawns the matching
# script of ../download and streams its output; GET /api/met/download/status polls it;
# POST /api/met/download/cancel kills it.  The scripts write the AVAILABLE index, so a
# run on the new directory starts without scanning.  They are Julia scripts run with the
# easlee project, so no Python is needed; the Earthdata token is passed through the
# environment, never on the command line.
const DOWNLOAD_DIR = normpath(joinpath(@__DIR__, "..", "download"))
const DL = Dict{String,Any}("running" => false, "log" => String[], "exit" => nothing,
                            "proc" => nothing, "cmd" => "", "out" => "", "model" => "")
const DL_LOCK = ReentrantLock()
const EASLEE_PROJECT = normpath(joinpath(@__DIR__, "..", "easlee"))
_julia() = joinpath(Sys.BINDIR, Sys.iswindows() ? "julia.exe" : "julia")
function _dl_log!(line::AbstractString)
    lock(DL_LOCK) do
        for piece in split(line, '\r')            # progress lines use \r
            isempty(strip(piece)) || push!(DL["log"], String(piece))
        end
        length(DL["log"]) > 3000 && deleteat!(DL["log"], 1:1000)
    end
end
_s(b, k, d="") = haskey(b, k) && b[k] !== nothing ? String(string(b[k])) : d
function _download_command(b)
    model = uppercase(_s(b, "model", "GFS")); out = _s(b, "out")
    isempty(out) && error("target directory is empty")
    env = copy(ENV)
    if model == "GFS"
        args = ["--out", out, "--res", _s(b, "res", "0p25"), "--source", _s(b, "source", "auto")]
        if _s(b, "mode", "forecast") == "forecast"
            isempty(_s(b, "analysis")) && error("analysis time is empty")
            append!(args, ["--analysis", _s(b, "analysis"), "--hours", _s(b, "hours", "48")])
        else
            isempty(_s(b, "start")) && error("start date is empty")
            append!(args, ["--start", _s(b, "start"), "--end", _s(b, "end", _s(b, "start")),
                           "--cycles", _s(b, "cycles", "00,12"), "--steps", _s(b, "steps", "3,6,9,12,15")])
        end
        script = "gfs.jl"
    elseif model == "MERRA2"
        isempty(_s(b, "start")) && error("start date is empty")
        args = ["--start", _s(b, "start"), "--end", _s(b, "end", _s(b, "start")), "--out", out]
        isempty(_s(b, "collections")) || append!(args, ["--collections", _s(b, "collections")])
        isempty(_s(b, "token")) || (env["EARTHDATA_TOKEN"] = _s(b, "token"))
        script = "merra2.jl"
    elseif model == "ERA5"
        args = ["--out", out]
        if !isempty(_s(b, "url"))
            append!(args, ["--url", _s(b, "url")])
        else
            append!(args, ["--dataset", lowercase(_s(b, "dataset", "etex"))])
        end
        script = "era5_archive.jl"
    else
        error("unknown model $model (GFS, MERRA2 or ERA5)")
    end
    cmd = Cmd([_julia(), "--project=$EASLEE_PROJECT", joinpath(DOWNLOAD_DIR, script), args...])
    return setenv(cmd, env), "$(script) $(join(args, ' '))", model, out
end
function handle_met_download(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    DL["running"] && return json_err("a download is already running")
    body = JSON3.read(String(req.body))
    local cmd, shown, model, out
    try
        cmd, shown, model, out = _download_command(body)
    catch err
        return json_err(sprint(showerror, err))
    end
    isdir(DOWNLOAD_DIR) || return json_err("download scripts not found: $DOWNLOAD_DIR")
    lock(DL_LOCK) do
        empty!(DL["log"]); DL["exit"] = nothing; DL["cmd"] = shown; DL["out"] = out; DL["model"] = model
        DL["running"] = true
    end
    _dl_log!("\$ " * shown)
    Threads.@spawn begin
        try
            proc = open(pipeline(cmd; stderr=stdout))
            DL["proc"] = proc
            for line in eachline(proc)
                _dl_log!(line)
            end
            wait(proc)
            DL["exit"] = proc.exitcode
            _dl_log!(proc.exitcode == 0 ? "finished — AVAILABLE index written in $out" :
                     "exited with code $(proc.exitcode)")
        catch err
            DL["exit"] = -1; _dl_log!("error: " * sprint(showerror, err))
        finally
            DL["running"] = false; DL["proc"] = nothing
        end
    end
    return json_ok(Dict("ok" => true, "cmd" => shown))
end
function handle_met_download_status(::HTTP.Request)
    lines = lock(DL_LOCK) do
        n = length(DL["log"]); DL["log"][max(1, n - 199):n]
    end
    json_ok(Dict("running" => DL["running"], "exit" => DL["exit"], "cmd" => DL["cmd"],
                 "out" => DL["out"], "model" => DL["model"], "log" => lines,
                 "julia" => _julia(), "scripts" => isdir(DOWNLOAD_DIR)))
end
function handle_met_download_cancel(req::HTTP.Request)
    g = _post_guard(req); g === nothing || return g
    p = DL["proc"]
    (p === nothing || !DL["running"]) && return json_err("no download running")
    try; kill(p); catch; end
    _dl_log!("cancelled")
    json_ok(Dict("ok" => true))
end

const ROUTER = HTTP.Router()
HTTP.register!(ROUTER, "GET",  "/",             handle_index)
HTTP.register!(ROUTER, "GET",  "/static/**",    handle_static)   # legacy
HTTP.register!(ROUTER, "GET",  "/app.js",       handle_asset)
HTTP.register!(ROUTER, "GET",  "/schema_form.js", handle_asset)
HTTP.register!(ROUTER, "GET",  "/style.css",    handle_asset)
HTTP.register!(ROUTER, "GET",  "/easlee_logo.svg", handle_asset)
HTTP.register!(ROUTER, "GET",  "/easlee_mark.svg", handle_asset)
HTTP.register!(ROUTER, "GET",  "/coastline.json", handle_asset)
HTTP.register!(ROUTER, "GET",  "/api/state",    handle_state)
handle_species(::HTTP.Request) = HTTP.Response(200, ["Content-Type"=>"application/json"],
    JSON3.write(Dict("species" => isempty(LIVE.last_nc) ? String[] : list_species_in_nc(LIVE.last_nc))))
HTTP.register!(ROUTER, "GET",  "/api/figure",   handle_figure)
HTTP.register!(ROUTER, "GET",  "/api/species",  handle_species)
HTTP.register!(ROUTER, "POST", "/api/control",  handle_control)
HTTP.register!(ROUTER, "POST", "/api/eval",     handle_eval)
HTTP.register!(ROUTER, "GET",  "/api/schema",     handle_schema)
HTTP.register!(ROUTER, "GET",  "/api/config",     handle_get_config)
HTTP.register!(ROUTER, "POST", "/api/config",     handle_set_config)
HTTP.register!(ROUTER, "POST", "/api/load_toml",  handle_load_toml)
HTTP.register!(ROUTER, "GET",  "/api/config.toml", handle_config_toml)
HTTP.register!(ROUTER, "GET",  "/api/browse",     handle_browse)
HTTP.register!(ROUTER, "GET",  "/api/domain",     handle_domain)
HTTP.register!(ROUTER, "POST", "/api/met/download",        handle_met_download)
HTTP.register!(ROUTER, "GET",  "/api/met/download/status", handle_met_download_status)
HTTP.register!(ROUTER, "POST", "/api/met/download/cancel", handle_met_download_cancel)

start_server(port::Int) = HTTP.serve(ROUTER, get(ENV, "EASLEE_GUI_HOST", "127.0.0.1"), port)   # loopback by default (review 2026-09-23)
