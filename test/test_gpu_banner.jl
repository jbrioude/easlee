# test_gpu_banner.jl — the startup GPU status line must say unambiguously
# whether THIS process uses the GPU, in every configuration:
#   (a) nothing requested            → "GPU: not requested — CPU run"
#   (b) requested, device unusable   → a warning naming the fallback
# (A functional-GPU line is exercised by hand — it needs hardware.)
# Motivated by a first-time HPC user (2026-09-04) who had the GPU allocated
# to the job but no Julia process on it: nvidia-smi showed the allocation,
# not the use.  Each case runs in its own process because the GPU gate is a
# load-time constant read from ENV.

const REPO = abspath(joinpath(@__DIR__, ".."))

npass = 0; nfail = 0
function check(name, ok; info = "")
    global npass, nfail
    if ok
        npass += 1; println("PASS  ", name)
    else
        nfail += 1; println("FAIL  ", name, isempty(info) ? "" : "  ($info)")
    end
end

loader = tempname() * ".jl"
write(loader, """
const REPO = "$(escape_string(REPO))"
global gpu = false; global chem = false
include(joinpath(REPO, "easlee", "src", "config", "EasleeConfig.jl")); using .EasleeConfig
include(joinpath(REPO, "easlee", "src", "io", "initialize.jl"))
include(joinpath(REPO, "easlee", "src", "io", "tools.jl"))
include(joinpath(REPO, "easlee", "src", "io", "output.jl"))
include(joinpath(REPO, "easlee", "EASLEE.jl")); using .EASLEE
EASLEE.Operators.gpu_status_banner()
""")

jl = `$(Base.julia_cmd()) --project=$(joinpath(REPO, "easlee")) $loader`

function run_case(env)
    out = IOBuffer()
    ok = withenv(env...) do
        success(run(pipeline(jl; stdout = out, stderr = out); wait = true))
    end
    ok, String(take!(out))
end

# (a) no request → CPU line, no warning
ok, txt = run_case(["EASLEE_ADVANCE" => nothing, "EASLEE_GPU" => nothing])
check("cpu: process exits 0", ok)
check("cpu: 'GPU: not requested' line printed", occursin("GPU: not requested", txt); info = txt)
check("cpu: says how to enable it", occursin("EASLEE_ADVANCE=ka_gpu", txt))
check("cpu: no GPU warning", !occursin("GPU requested", txt))

# (b) requested but no usable device (CUDA hidden, or not installed at all)
ok, txt = run_case(["EASLEE_ADVANCE" => "ka_gpu", "EASLEE_GPU" => nothing,
                    "CUDA_VISIBLE_DEVICES" => ""])
check("nodev: process exits 0 (fallback, never a crash)", ok)
check("nodev: warning names the request", occursin("GPU requested (EASLEE_ADVANCE=ka_gpu)", txt); info = txt)
check("nodev: warning states the CPU fallback", occursin("CPU run", txt))
check("nodev: no false 'GPU backend enabled'", !occursin("GPU backend enabled", txt))

rm(loader; force = true)
println("\n", npass, " passed, ", nfail, " failed")
nfail == 0 || exit(1)
