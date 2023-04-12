module ComplexPaths

export AbstractPath, Path, CirclePath, pointsonpath

"""
    Path(parameterization::Function, start::Complex, end::Complex)

    # Feilds
    - `parameterization::Function`: A parametric function of a complex path.
    - `start::Complex`: Initial value for the parametric function.
    - `end::Complex`: Final value for the parametric function.
"""

abstract type AbstractPath end

struct Path <: AbstractPath
    parameterization::Function
    start::Complex
    ending::Complex
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

end # module ComplexPaths

