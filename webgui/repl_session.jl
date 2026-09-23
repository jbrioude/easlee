# Simple eval-loop REPL: captures the result repr; print/println output goes to
# the terminal running Julia (stdout capture can be added later if needed).
function eval_expression(code::String) :: Tuple{String, Bool}
    code = strip(code)
    isempty(code) && return ("", false)

    # Multi-expression: handle `1+1; 2+2` by wrapping in a begin block
    src = occursin('\n', code) ? "begin\n$code\nend" : code

    try
        expr   = Meta.parse(src)
        result = Core.eval(Main, expr)
        result === nothing && return ("", false)
        return (repr(result), false)
    catch e
        return (sprint(showerror, e), true)
    end
end
