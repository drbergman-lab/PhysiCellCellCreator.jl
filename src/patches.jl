# ---------------------------------------------------------------------------
# Patch type registry
# ---------------------------------------------------------------------------

struct DiscPatch end
struct AnnulusPatch end
struct RectanglePatch end
struct EverywherePatch end
struct EllipticalDiscPatch end
struct EllipsePatch end

const patch_type_registry = [
    "disc" => DiscPatch,
    "annulus" => AnnulusPatch,
    "rectangle" => RectanglePatch,
    "everywhere" => EverywherePatch,
    "elliptical_disc" => EllipticalDiscPatch,
    "ellipse" => EllipsePatch,
]

"""
    supportedPatchTypes()

The values accepted by the `type` attribute of a `patch_collection` element.
"""
supportedPatchTypes() = first.(patch_type_registry)

function patchType(s::AbstractString)
    ind = findfirst(p -> p.first == s, patch_type_registry)
    @assert !isnothing(ind) "Patch type $(s) is not supported. Supported types are: $(supportedPatchTypes())"
    return patch_type_registry[ind].second
end

const centered_frame_elements = ["normal", "rotation"]
const placement_elements = ["number", "spacing", "max_fails", "restrict_to_domain"]

"""
    geometryElements(T)

The elements describing the size, position, and orientation of a patch of type `T`; these are also the elements of the carveout of the same type.
"""
geometryElements(::Type{DiscPatch}) = ["x0", "y0", "z0", "radius", centered_frame_elements...]
geometryElements(::Type{AnnulusPatch}) = ["x0", "y0", "z0", "inner_radius", "outer_radius", centered_frame_elements...]
geometryElements(::Type{EllipticalDiscPatch}) = ["x0", "y0", "z0", "a", "b", centered_frame_elements...]
geometryElements(::Type{EllipsePatch}) = ["x0", "y0", "z0", "a", "b", "thickness", centered_frame_elements...]
geometryElements(::Type{RectanglePatch}) = ["x0", "y0", "z0", "x1", "y1", "width", "height"]

"""
    supportedPatchTextElements(patch_type)

The elements a `patch` of the given type may contain.
"""
supportedPatchTextElements(::Type{EverywherePatch}) = ["number", "spacing", "max_fails"]
supportedPatchTextElements(::Type{T}) where {T} = [geometryElements(T); placement_elements]
supportedPatchTextElements(patch_type::AbstractString) = supportedPatchTextElements(patchType(patch_type))

# ---------------------------------------------------------------------------
# Geometry: XML -> (Shape, Frame)
# ---------------------------------------------------------------------------

parseCenteredFrame(patch::XMLElement, descriptor::String) = Frame(parseCenter(patch, descriptor), parseNormal(patch, descriptor), parseRotation(patch, descriptor))

"""
    parseGeometry(T, patch, descriptor)

Parse the geometric elements of a patch (or carveout) of type `T` into a [`Shape`](@ref) and its [`Frame`](@ref).
"""
function parseGeometry(::Type{DiscPatch}, patch::XMLElement, descriptor::String)
    radius = parsePositiveFloatElement(patch, "radius", descriptor)
    return EllipticalDisc(radius, radius), parseCenteredFrame(patch, descriptor)
end

function parseGeometry(::Type{AnnulusPatch}, patch::XMLElement, descriptor::String)
    inner_radius = parseFloatElement(patch, "inner_radius", descriptor)
    outer_radius = parsePositiveFloatElement(patch, "outer_radius", descriptor)
    @assert inner_radius >= 0 "inner_radius for $(descriptor) must be non-negative. Got $(inner_radius)."
    @assert inner_radius <= outer_radius "Inner radius must be less than or equal to outer radius. Got inner_radius=$(inner_radius) > outer_radius=$(outer_radius)."
    shape = inner_radius == 0 ? EllipticalDisc(outer_radius, outer_radius) : EllipticalRing(inner_radius, inner_radius, outer_radius - inner_radius)
    return shape, parseCenteredFrame(patch, descriptor)
end

function parseGeometry(::Type{EllipticalDiscPatch}, patch::XMLElement, descriptor::String)
    a = parsePositiveFloatElement(patch, "a", descriptor)
    b = parsePositiveFloatElement(patch, "b", descriptor)
    return EllipticalDisc(a, b), parseCenteredFrame(patch, descriptor)
end

function parseGeometry(::Type{EllipsePatch}, patch::XMLElement, descriptor::String)
    a = parsePositiveFloatElement(patch, "a", descriptor)
    b = parsePositiveFloatElement(patch, "b", descriptor)
    thickness = parsePositiveFloatElement(patch, "thickness", descriptor)
    return EllipticalRing(a, b, thickness), parseCenteredFrame(patch, descriptor)
end

function parseGeometry(::Type{RectanglePatch}, patch::XMLElement, descriptor::String)
    center, width, height = parseRectangleParameters(patch, descriptor)
    return Rectangle(width, height), Frame(center, [0.0, 0.0, 1.0], 0.0)
end

# ---------------------------------------------------------------------------
# Placement settings
# ---------------------------------------------------------------------------

"""
    CellRadius

Accepted types for the `cell_radius` keyword argument of [`generateICCell`](@ref): `nothing`, a single `Real`, or a `Dict` mapping cell type names to radii.
"""
const CellRadius = Union{Nothing,Real,AbstractDict}

"""
    HexagonalDefaults

Caller-supplied settings for hexagonal grid patches: the cell radius (per cell type or global) and the spacing factor, giving the default `spacing = spacing_factor * 2 * cell_radius` for patches that omit a `spacing` element, and the boundary tolerance (as a fraction of the spacing) used when deciding whether a lattice point lies on the boundary of a patch.
"""
struct HexagonalDefaults
    cell_radius::CellRadius
    spacing_factor::Float64
    boundary_tolerance::Float64
end

resolveCellRadius(::Nothing, ::String) = nothing
resolveCellRadius(cell_radius::Real, ::String) = Float64(cell_radius)
function resolveCellRadius(cell_radius::AbstractDict, cell_type::String)
    haskey(cell_radius, cell_type) || return nothing
    radius = cell_radius[cell_type]
    @assert radius isa Real "cell_radius for $(cell_type) must be a Real. Got $(radius) of type $(typeof(radius))."
    return Float64(radius)
end

"""
    supportedPlacements()

The values accepted by the `placement` attribute of a `patch` element.
`random` (the default when the attribute is absent) places `number` cells uniformly at random; `hexagonal` fills the patch with a hexagonal grid of cells with the given `spacing`.
"""
supportedPlacements() = ["random", "hexagonal"]

function parsePlacement(patch::XMLElement, descriptor::String)
    placement = attribute(patch, "placement")
    if isnothing(placement)
        return "random"
    end
    @assert placement in supportedPlacements() "placement=\"$(placement)\" for $(descriptor) is not supported. Supported placements are: $(supportedPlacements())"
    return placement
end

"""
    parseHexagonalGrid(patch, cell_type, hex_defaults)

Return `nothing` if the patch's `placement` attribute is not `hexagonal`; otherwise, return the lattice spacing (distance between neighboring cell centers).

The spacing comes from the patch's `spacing` element if present, and otherwise from `hex_defaults.spacing_factor * 2 * cell_radius`.
A hexagonal patch cannot have a `number` element since the patch size and spacing determine how many cells are placed, and a random patch cannot have a `spacing` element.
"""
function parseHexagonalGrid(patch::XMLElement, cell_type::String, hex_defaults::HexagonalDefaults)
    descriptor = patchDescriptor(patch, cell_type)
    placement = parsePlacement(patch, descriptor)
    spacing_element = find_element(patch, "spacing")
    if placement == "random"
        @assert isnothing(spacing_element) "spacing is set for $(descriptor) but its placement is random. spacing is only used with placement=\"hexagonal\"."
        return nothing
    end
    @assert isnothing(find_element(patch, "number")) "number is set for $(descriptor) but its placement is hexagonal. These are incompatible as both determine the number of cells placed in the patch."
    if !isnothing(spacing_element)
        spacing = parse(Float64, strip(content(spacing_element)))
    else
        radius = resolveCellRadius(hex_defaults.cell_radius, cell_type)
        @assert !isnothing(radius) "placement is hexagonal for $(descriptor) but no spacing element is set and no cell_radius was passed to generateICCell for this cell type."
        spacing = hex_defaults.spacing_factor * 2 * radius
    end
    @assert spacing > 0 "spacing for $(descriptor) must be positive. Got $(spacing)."
    return spacing
end

"""
    max_hexagonal_grid_points

Upper bound on the number of lattice points a single hexagonal grid patch may generate before `hexagonalGridPoints` throws, guarding against a mistyped `spacing` exhausting memory.
"""
const max_hexagonal_grid_points = 10_000_000

"""
    hexagonalGridPoints(spacing, half_width, half_height, boundary_tolerance)

Return the `(c1, c2)` coordinates of a hexagonal lattice with the given `spacing` covering `[-half_width, half_width] x [-half_height, half_height]`.

The lattice has a point at the origin, rows parallel to the first axis separated by `spacing * sqrt(3) / 2`, and alternating rows offset by `spacing / 2`.
Points within `boundary_tolerance * spacing` outside the region are still kept so floating-point rounding does not drop boundary points.
An `ArgumentError` is thrown if the lattice would exceed `max_hexagonal_grid_points` points, which almost certainly indicates a spacing typo.
"""
function hexagonalGridPoints(spacing::Float64, half_width::Float64, half_height::Float64, boundary_tolerance::Float64)
    row_height = spacing * sqrt(3) / 2
    n_rows = ceil(Int, half_height / row_height)
    n_cols = ceil(Int, half_width / spacing) + 1
    n_points_bound = (2 * n_rows + 1) * (2 * n_cols + 1)
    if n_points_bound > max_hexagonal_grid_points
        throw(ArgumentError("A hexagonal grid with spacing $(spacing) over a $(2 * half_width) x $(2 * half_height) region would have up to $(n_points_bound) points, exceeding the limit of $(max_hexagonal_grid_points). Check the spacing."))
    end
    tol = boundary_tolerance * spacing
    c1 = Float64[]
    c2 = Float64[]
    for j in -n_rows:n_rows
        y = j * row_height
        abs(y) <= half_height + tol || continue
        offset = isodd(j) ? spacing / 2 : 0.0
        for i in -n_cols:n_cols
            x = i * spacing + offset
            abs(x) <= half_width + tol || continue
            push!(c1, x)
            push!(c2, y)
        end
    end
    return c1, c2
end

# ---------------------------------------------------------------------------
# Placing cells in a shape
# ---------------------------------------------------------------------------

function createCellsDataFrame(number::Int, cell_coords_fn::Function, restrict_to_domain::Bool, domain_dict::Dict{<:AbstractString,<:Real}, carveouts::Vector{Carveout}, max_fails::Int)
    total_placed = 0
    df = DataFrame(x=Float64[], y=Float64[], z=Float64[])
    fails_remaining = max_fails
    while total_placed < number
        df_new = cell_coords_fn(number - total_placed)
        filterCells!(df_new, restrict_to_domain, domain_dict, carveouts)
        append!(df, df_new)
        total_placed += size(df_new, 1)
        if size(df_new, 1) == 0
            if fails_remaining == 0
                throw(ArgumentError("Failed to place a single cell after $(max_fails) consecutive attempts."))
            end
            fails_remaining -= 1
        else
            fails_remaining = max_fails
        end
    end
    return df
end

function warnIfEmpty(cell_df::DataFrame, descriptor::String)
    if isempty(cell_df)
        @warn "Hexagonal grid for $(descriptor) placed no cells. Check the spacing, patch size, domain, and carveouts."
    end
end

"""
    placeShape!(df, shape, frame, restrict_to_domain, patch, cell_type, domain_dict, hex_defaults)

Place cells of `cell_type` in `shape` (positioned by `frame`), either uniformly at random or on a hexagonal grid centered on the shape, then append them to `df`.
"""
function placeShape!(df::DataFrame, shape::PlanarShape, frame::Frame, restrict_to_domain::Bool, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    descriptor = patchDescriptor(patch, cell_type)
    carveouts = parseCarveouts(patch, descriptor)
    spacing = parseHexagonalGrid(patch, cell_type, hex_defaults)

    if isnothing(spacing)
        number = parseNumber(patch, descriptor)
        max_fails = parseMaxFails(patch, descriptor)
        cell_coords_fn(n::Int) = toDomain(frame, sampleUniform(shape, n)...)
        cell_df = createCellsDataFrame(number, cell_coords_fn, restrict_to_domain, domain_dict, carveouts, max_fails)
    else
        half_width, half_height = bounds(shape)
        c1, c2 = hexagonalGridPoints(spacing, half_width, half_height, hex_defaults.boundary_tolerance)
        keep_ind = inside(shape, c1, c2; tol=hex_defaults.boundary_tolerance * spacing)
        cell_df = toDomain(frame, c1[keep_ind], c2[keep_ind])
        filterCells!(cell_df, restrict_to_domain, domain_dict, carveouts)
        warnIfEmpty(cell_df, descriptor)
    end
    cell_df[!, :type] .= cell_type
    append!(df, cell_df)
end

function generatePatch!(::Type{T}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults) where {T}
    descriptor = patchDescriptor(patch, cell_type)
    shape, frame = parseGeometry(T, patch, descriptor)
    restrict_to_domain = parseRestrictToDomain(patch, descriptor)
    placeShape!(df, shape, frame, restrict_to_domain, patch, cell_type, domain_dict, hex_defaults)
end

function generatePatch!(::Type{EverywherePatch}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    # everywhere patches are 2D: cells are placed at z = 0 regardless of the domain z bounds (both random and hexagonal placement)
    width = Float64(domain_dict["x_max"] - domain_dict["x_min"])
    height = Float64(domain_dict["y_max"] - domain_dict["y_min"])
    center = [domain_dict["x_min"] + width / 2, domain_dict["y_min"] + height / 2, 0.0]
    frame = Frame(center, [0.0, 0.0, 1.0], 0.0)
    placeShape!(df, Rectangle(width, height), frame, false, patch, cell_type, domain_dict, hex_defaults)
end
