module ComplexPaths

using Intervals

import Intervals: (..)

export AbstractPath, Path, ClosedPath, CirclePath, PiecewisePath, pointsonpath, (..)

abstract type AbstractPath end

"""@docs
Path(parameterization::Function, domain::AbstractInterval)

# Feilds
- `parameterization::Function`: A parametric function of a complex path.
- `start::Real`: Initial value for the parametric function.
- `end::Real`: Final value for the parametric function.
"""
struct Path <: AbstractPath
    parameterization::Function
    domain::AbstractInterval
end

"""@docs
ClosedPath(parameterization::Function, domain::AbstractInterval)

Same as Path with added closure assertion.
"""
struct ClosedPath <: AbstractPath
    parameterization::Function
    domain::AbstractInterval
    function ClosedPath(parameterization,domain) 
        path = Path(parameterization,domain)
        @assert isclosed(path)
        return new(parameterization,domain)
    end
end

#Map values linaerly from X:[a,b] -> Y:[c,d]
function _map_values(t,X::Interval,Y::Interval)::Real

    a = X.first
    b = X.last
    c = Y.first 
    d = Y.last

    t′ = c + ((d-c)/(b-a))*(t-a)

    return t′

end

function _intervaldist(α::AbstractPath)
    D = α.domain
    lesser, greater = D.first < D.last ? (D.first, D.last) : (D.last, D.first)
    return abs(greater - lesser)
end

#Build a dictionary mapping a range spanning all paths into their constituant domains.
function _build_piecewise_dict(P1, Pn...)
    
    out = Dict{Interval,Path}()

    P_mag = _intervaldist(P1)
    out[0.0..P_mag] = P1
    lv = P_mag
    
    if length(Pn) ≥ 1
        for P ∈ Pn
            P_mag = _intervaldist(P)
            pair = Interval{Open,Closed}(lv,lv+P_mag) => P
            push!(out,pair)
           
            lv += P_mag
        end
    end

    return lv, out

end
  
function _map_call(t,key,value)
    t_map = _map_values(t, key, value.domain)
    f = value.parameterization
    return f(t_map)
end

function _call_piecewise(t::Real, D::AbstractDict)
    for (key,value) ∈ D
        if t ∈ key
            return key, _map_call(t,key,value)
        end
    end
end

function _call_piecewise(t::Real, D::AbstractDict, previous::Interval)
    if t ∈ previous
        value = D[previous]
        return previous, _map_call(t,previous,value)
    else
        return _call_piecewise(t,D)
    end
end

"""@docs
PiecewisePath(parameterization::Function, domain::AbstractInterval)

A path composed of several parameterizations mapped across one domain.
"""
struct PiecewisePath <: AbstractPath
    parameterization::Function
    domain::AbstractInterval
    piecewise::Dict{Interval,Path}
    _call_piecewise_internal(t) = _call_piecewise(t,piecewise)
    _call_piecewise_internal(t,prev) = _call_piecewise(t,piecewise,prev)
    function PiecewisePath(Ps...)
        domain_size, piecewise = _build_piecewise_dict(Ps...)
        return new(_call_piecewise_internal,0.0..domain_size,piecewise)
    end
    function PiecewisePath(parameterization,domain)
        domain_size, piecewise = _build_piecewise_dict(Path(parameterization,domain))
        return new(_call_piecewise_internal,0.0..domain_size,piecewise)
    end
end

"""@docs
CirclePath(r::Real)

Define a circular path of radius `r` centered around the origin. 
"""
CirclePath(r::Real)::Path = Path(t -> (r*ℯ^(t*1im)),0..2π)

"""@docs
pointsonpath(P::AbstractPath, n::Integer)

Return a vector of `n` evenly spaced points along path `P`.
"""
function pointsonpath(P::AbstractPath, n::Integer)::Vector
    @assert n > 0 "n must be a positive non-zero integer"
    D = P.domain
    return [P.parameterization(t) for t in range(D.first, D.last, length=n)]
end 

function pointsonpath(P::PiecewisePath, n::Integer)::Vector

    @assert n > 0 "n must be a positive non-zero integer"
    D = P.domain
    prev = (first ∘ keys)(P.piecewise)
    
    out = []

    for t ∈ range(D.first, D.last, length=n)
        prev, val = P.parameterization(t,prev)
        push!(out,val)
    end

    return out
end

"""@docs
*(α::AbstractPath, β::AbstractPath)::PiecewisePath

Concatenate two paths together, yeilding a PiecewisePath
"""
function Base.:*(α::AbstractPath, β::AbstractPath)::PiecewisePath

    return PiecewisePath(α,β)
end

function Base.:*(P1::AbstractPath, Pn::AbstractPath...)::PiecewisePath

    return PiecewisePath(P1, Pn...)
end

"""@docs
isclosed(γ::AbstractPath)::Bool

Check if `γ` is a closed Path. i.e. the start and end points are the same (within a tolerance because of floating point math).
"""
function isclosed(γ::AbstractPath)::Bool
    return γ.parameterization(γ.start) ≈ γ.parameterization(γ.ending)
end

end # module ComplexPaths

