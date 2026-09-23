
# Vertical level types.
# Each type describes a different vertical coordinate system;
# dispatch on ztoindl() selects the right conversion.

abstract type AbstractLevel end

struct SigmaPressureLevel <: AbstractLevel   # AROME, ECMWF
    Nk::Int32
    ak::Vector{Float32}
    bk::Vector{Float32}
end

mutable struct SigmaGeopotentialLevel <: AbstractLevel  # MesoNH
    Nk::Int32
    wname::String
    zname::String
    Zlevels::Vector{Float32}
end

mutable struct PressureLevel <: AbstractLevel            # ECMWF, GFS
    Nk::Int32
    wname::String
    zname::String
    plevels::Vector{Float32}
end

mutable struct HHLLevel <: AbstractLevel                 # ICON
    Nk::Int32
    wname::String
    zname::String
    hhl::Array{Float32,3}
end

mutable struct PotentialTemperatureLevel <: AbstractLevel  # diabatic
    Nk::Int32
    wname::String
    zname::String
    thetalevels::Vector{Float32}
end
