"""
    Shape

A region in which cells are placed (a patch) or from which they are removed (a carveout), defined in the local coordinates of a [`Frame`](@ref).

Every shape implements [`sampleUniform`](@ref), [`inside`](@ref), and [`bounds`](@ref).
"""
abstract type Shape end

"""
    PlanarShape

A two-dimensional [`Shape`](@ref) lying in the plane `c3 = 0` of its frame.
Its methods take and return the in-plane coordinates `(c1, c2)` only.
"""
abstract type PlanarShape <: Shape end

"""
    EllipticalDisc(a, b)

The filled ellipse `(c1 / a)^2 + (c2 / b)^2 <= 1`. A disc is the special case `a == b`.
"""
struct EllipticalDisc <: PlanarShape
    a::Float64
    b::Float64
    function EllipticalDisc(a::Real, b::Real)
        @assert a > 0 && b > 0 "The axes of an elliptical disc must be positive. Got a=$(a), b=$(b)."
        return new(a, b)
    end
end

"""
    EllipticalRing(a, b, thickness)

The region between the ellipse with axes `(a, b)` and the ellipse with axes `(a + thickness, b + thickness)`, following the `ellipse` patch of PhysiCellECMCreator.jl.
An annulus is the special case `a == b`.
"""
struct EllipticalRing <: PlanarShape
    a::Float64
    b::Float64
    thickness::Float64
    function EllipticalRing(a::Real, b::Real, thickness::Real)
        @assert a > 0 && b > 0 "The axes of an ellipse must be positive. Got a=$(a), b=$(b)."
        @assert thickness >= 0 "The thickness of an ellipse must be non-negative. Got $(thickness)."
        @assert thickness > 0 || a == b "The thickness of a non-circular ellipse must be positive (a zero-thickness ellipse has no area). Got a=$(a), b=$(b), thickness=$(thickness)."
        return new(a, b, thickness)
    end
end

"""
    Rectangle(width, height)

The axis-aligned rectangle `[-width/2, width/2] x [-height/2, height/2]` centered at the origin of its frame.
"""
struct Rectangle <: PlanarShape
    width::Float64
    height::Float64
    function Rectangle(width::Real, height::Real)
        @assert width > 0 && height > 0 "The width and height of a rectangle must be positive. Got width=$(width), height=$(height)."
        return new(width, height)
    end
end

ellipseFunctional(c1::Real, c2::Real, a::Real, b::Real) = (c1 / a)^2 + (c2 / b)^2

"""
    max_proposal_batch

Upper bound on the number of candidate points proposed in one round of rejection sampling, bounding temporary memory for very thin rings (which then take more rounds).
"""
const max_proposal_batch = 1_000_000

"""
    inside(shape, c1, c2; tol=0.0)

Boolean vector saying which local points `(c1, c2)` lie in `shape`, including its boundary.
This is the single membership test used by carveouts, lattice filtering, and rejection sampling.

`tol` enlarges the shape (for a ring, it also shrinks the inner ellipse) so that lattice points which land on the boundary up to floating-point rounding are kept.
It is only passed for hexagonal placement; carveouts and random sampling use `tol = 0`, since random points hit a boundary with probability zero.
For an ellipse, `tol` is added to the axes, which keeps a slightly thinner band than a true distance band would; that is immaterial at the rounding scale it is meant for.
"""
function inside(shape::EllipticalDisc, c1::AbstractVector{<:Real}, c2::AbstractVector{<:Real}; tol::Real=0.0)
    return ellipseFunctional.(c1, c2, shape.a + tol, shape.b + tol) .<= 1
end

function inside(shape::EllipticalRing, c1::AbstractVector{<:Real}, c2::AbstractVector{<:Real}; tol::Real=0.0)
    a_outer = shape.a + shape.thickness + tol
    b_outer = shape.b + shape.thickness + tol
    in_outer = ellipseFunctional.(c1, c2, a_outer, b_outer) .<= 1
    a_inner = shape.a - tol
    b_inner = shape.b - tol
    if a_inner <= 0 || b_inner <= 0
        return in_outer
    end
    strictly_in_inner = ellipseFunctional.(c1, c2, a_inner, b_inner) .< 1
    return in_outer .& .!strictly_in_inner
end

function inside(shape::Rectangle, c1::AbstractVector{<:Real}, c2::AbstractVector{<:Real}; tol::Real=0.0)
    return (abs.(c1) .<= shape.width / 2 + tol) .& (abs.(c2) .<= shape.height / 2 + tol)
end

"""
    bounds(shape)

Half-width and half-height of the smallest axis-aligned box (in local coordinates) containing `shape`, used to generate a lattice over the shape.
"""
bounds(shape::EllipticalDisc) = (shape.a, shape.b)
bounds(shape::EllipticalRing) = (shape.a + shape.thickness, shape.b + shape.thickness)
bounds(shape::Rectangle) = (shape.width / 2, shape.height / 2)

"""
    sampleUniform(shape, n)

Return the local coordinates `(c1, c2)` of `n` points sampled uniformly at random from `shape`.
"""
function sampleUniform(shape::EllipticalDisc, n::Int)
    r = sqrt.(rand(n))
    θ = 2π .* rand(n)
    return shape.a .* r .* cos.(θ), shape.b .* r .* sin.(θ)
end

function sampleUniform(shape::EllipticalRing, n::Int)
    a, b, t = shape.a, shape.b, shape.thickness
    if a == b
        # circular annulus: sample the radius exactly
        r = sqrt.(a^2 .+ ((a + t)^2 - a^2) .* rand(n))
        θ = 2π .* rand(n)
        return r .* cos.(θ), r .* sin.(θ)
    end
    # the outer ellipse is an offset of the inner one, not a scaling, so sample the outer disc and reject the inside of the inner ellipse
    outer = EllipticalDisc(a + t, b + t)
    acceptance = 1 - a * b / ((a + t) * (b + t))
    c1 = Float64[]
    c2 = Float64[]
    while length(c1) < n
        # propose enough candidates that one round usually suffices: 20% more than the expected need absorbs the
        # binomial fluctuation in how many are accepted, and the +8 does the same for small n, where that fluctuation
        # is large relative to n. Neither value is special; a second round just proposes the shortfall.
        m = min(ceil(Int, 1.2 * (n - length(c1)) / acceptance) + 8, max_proposal_batch)
        d1, d2 = sampleUniform(outer, m)
        # no tolerance: it only matters for lattice points landing exactly on the boundary, which random samples do with probability zero
        keep_ind = inside(shape, d1, d2)
        append!(c1, d1[keep_ind])
        append!(c2, d2[keep_ind])
    end
    return c1[1:n], c2[1:n]
end

function sampleUniform(shape::Rectangle, n::Int)
    return shape.width .* (rand(n) .- 0.5), shape.height .* (rand(n) .- 0.5)
end
