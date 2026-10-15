"""
    Frame(center, normal, rotation)

Local coordinate frame of a patch or carveout: an origin at `center` and an orthonormal basis `(u1, u2, n)` stored as the columns of `R`.

`n` is the unit `normal`.
`u1` is the unit projection of the x-axis onto the plane perpendicular to `n` (or of the y-axis if `n` is parallel to the x-axis), and `u2 = n × u1`.
Both are then rotated counter-clockwise about `n` by `rotation` radians.
For the default normal `(0, 0, 1)` and zero rotation, the frame is the domain's own axes.
The frame varies continuously with the normal except at normals parallel to the x-axis, where the fallback applies.

Shapes are defined in local coordinates `(c1, c2, c3)` along `(u1, u2, n)`; see [`toDomain`](@ref) and [`toLocal`](@ref).
"""
struct Frame
    center::Vector{Float64}
    R::Matrix{Float64}
end

function Frame(center::AbstractVector{<:Real}, normal::AbstractVector{<:Real}, rotation::Real)
    @assert length(center) == 3 "A frame center must have three coordinates. Got $(center)."
    @assert length(normal) == 3 "A normal vector must have three components. Got $(normal)."
    n = Float64.(normal)
    magnitude = norm(n)
    @assert magnitude > 0 "A normal vector must be non-zero."
    n ./= magnitude
    # n × (e × n) is the projection of e onto the plane perpendicular to n, computed without cancellation for nearly parallel vectors
    u1 = cross(n, cross([1.0, 0.0, 0.0], n))
    if norm(u1) < 1e-12 # the normal is parallel to the x-axis
        u1 = cross(n, cross([0.0, 1.0, 0.0], n))
    end
    u1 ./= norm(u1)
    u2 = cross(n, u1)
    c, s = cos(rotation), sin(rotation)
    v1 = c .* u1 .+ s .* u2
    v2 = -s .* u1 .+ c .* u2
    return Frame(Float64.(center), [v1 v2 n])
end

"""
    frameAlong(center, direction)

A [`Frame`](@ref) whose local x direction is `direction`, for shapes (such as a [`Segment`](@ref)) that only use their local x axis.
The normal is an arbitrary unit vector perpendicular to `direction` (perpendicular to the z-axis when possible), and the local y axis completes the right-handed system.
"""
function frameAlong(center::AbstractVector{<:Real}, direction::AbstractVector{<:Real})
    @assert length(center) == 3 "A frame center must have three coordinates. Got $(center)."
    @assert length(direction) == 3 "A direction vector must have three components. Got $(direction)."
    u1 = Float64.(direction)
    magnitude = norm(u1)
    @assert magnitude > 0 "A direction vector must be non-zero."
    u1 ./= magnitude
    n = cross(u1, [0.0, 0.0, 1.0])
    if norm(n) < 1e-12 # the direction is parallel to the z-axis
        n = cross(u1, [1.0, 0.0, 0.0])
    end
    n ./= norm(n)
    u2 = cross(n, u1)
    return Frame(Float64.(center), [u1 u2 n])
end

"""
    toDomain(frame, c1, c2)

Map in-plane local coordinates `(c1, c2)` (vectors) of a planar shape to a `DataFrame` with the domain coordinates `x`, `y`, `z`.
"""
function toDomain(frame::Frame, c1::AbstractVector{<:Real}, c2::AbstractVector{<:Real})
    R = frame.R
    x = frame.center[1] .+ R[1, 1] .* c1 .+ R[1, 2] .* c2
    y = frame.center[2] .+ R[2, 1] .* c1 .+ R[2, 2] .* c2
    z = frame.center[3] .+ R[3, 1] .* c1 .+ R[3, 2] .* c2
    return DataFrame(x=x, y=y, z=z)
end

"""
    toLocal(frame, df)

Map the domain coordinates `x`, `y`, `z` of the rows of `df` to local coordinates, returned as the three vectors `(c1, c2, c3)`.
"""
function toLocal(frame::Frame, df::AbstractDataFrame)
    R = frame.R
    dx = df.x .- frame.center[1]
    dy = df.y .- frame.center[2]
    dz = df.z .- frame.center[3]
    c1 = R[1, 1] .* dx .+ R[2, 1] .* dy .+ R[3, 1] .* dz
    c2 = R[1, 2] .* dx .+ R[2, 2] .* dy .+ R[3, 2] .* dz
    c3 = R[1, 3] .* dx .+ R[2, 3] .* dy .+ R[3, 3] .* dz
    return c1, c2, c3
end
