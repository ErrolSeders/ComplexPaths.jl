module ComplexPaths

export AbstractPath, Path, CirclePath, pointsonpath

"""
    Path(parameterization::Function, start::Complex, end::Complex)

    # Feilds
    - `parameterization::Function`: A parametric function of a complex path.
    - `start::Number`: Initial value for the parametric function.
    - `end::Number`: Final value for the parametric function.
"""

abstract type AbstractPath end

struct Path <: AbstractPath
    parameterization::Function
    start::Number
    ending::Number
end

"""
    CirclePath(r::Real)

    Define a circular path of radius `r` centered around the origin. 
"""
CirclePath(r::Real)::Path = Path(t -> (r*ℯ^(t*1im)),0,2π)

"""
    pointsonpath(P::Path, n::Integer)

    Return a vector of `n` evenly spaced points along path `P`.
"""
function pointsonpath(P::AbstractPath, n::Integer)::Vector
    @assert n > 0 "n must be a positive non-zero integer"
    return [P.parameterization(t) for t in range(P.start, P.ending, length=n)]
end 

"""
    *(α::AbstractPath, β::AbstractPath)

    Concatenate two paths together
"""
function *(α::AbstractPath, β::AbstractPath)
    @assert α.ending == β.start "Paths must meet at endpoint."
    α_par = α.parameterization
    β_par = β.parameterization
    new_par = t -> t ≤ α.ending ? α.parameterization(t) : β.parameterization(t)
    return Path(new_par, α.start, β.ending)
end

"""
    isclosed(γ::AbstractPath)::Bool

    Check if `γ` is a closed Path. i.e. the start and end points are the same (within a tolerance for floating point nonsense).

"""
function isclosed(γ::AbstractPath)::Bool
    return γ.parameterization(γ.start) ≈ γ.parameterization(γ.ending) atol=0.0001
end

end # module ComplexPaths

