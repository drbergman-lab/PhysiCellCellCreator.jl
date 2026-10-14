module PhysiCellCellCreator

using LightXML, CSV, DataFrames, LinearAlgebra, Compat

export generateICCell, createICCellXMLTemplate

@compat public supportedPatchTypes, supportedPatchTextElements, supportedCarveoutTypes, supportedCarveoutTextElements

"""
    CellRadius

Accepted types for the `cell_radius` keyword argument of [`generateICCell`](@ref): `nothing`, a single `Real`, or a `Dict` mapping cell type names to radii.
"""
const CellRadius = Union{Nothing,Real,AbstractDict{<:AbstractString,<:Real}}

"""
    HexagonalDefaults

Caller-supplied defaults used by hexagonal grid patches that omit a `spacing` element: the cell radius (per cell type or global) and the spacing factor, giving `spacing = spacing_factor * 2 * cell_radius`.
"""
struct HexagonalDefaults
    cell_radius::CellRadius
    spacing_factor::Float64
end

resolveCellRadius(::Nothing, ::String) = nothing
resolveCellRadius(cell_radius::Real, ::String) = Float64(cell_radius)
resolveCellRadius(cell_radius::AbstractDict{<:AbstractString,<:Real}, cell_type::String) = haskey(cell_radius, cell_type) ? Float64(cell_radius[cell_type]) : nothing

"""
    generateICCell(path_to_ic_cell_xml, domain_dict; output=nothing, cell_radius=nothing, spacing_factor=0.95)

Generate initial cell positions from the IC cell XML file at `path_to_ic_cell_xml`.

`domain_dict` must contain the keys `"x_min"`, `"x_max"`, `"y_min"`, `"y_max"`, `"z_min"`, and `"z_max"`.

# Keyword arguments
- `output`: if `nothing` (default), return a `DataFrame` with columns `x`, `y`, `z`, and `type`; otherwise, write a CSV file to this path.
- `cell_radius`: the radius of the cells being placed, used to compute the default spacing for patches that use a hexagonal grid but do not specify a `spacing` element.
  Either a single `Real` applied to every cell type or a `Dict` mapping cell type names (the `name` attribute of `cell_patches`) to radii.
  Only required when a hexagonal grid patch omits `spacing`.
- `spacing_factor`: multiplies the cell diameter to get the default spacing, i.e., `spacing = spacing_factor * 2 * cell_radius`.
  The default of `0.95` matches the slightly compressed packing used in PhysiCell's sample projects; `1.0` places cells tangent to their neighbors.
"""
function generateICCell(path_to_ic_cell_xml::AbstractString, domain_dict::Dict{<:AbstractString,<:Real}; output::Union{Nothing,AbstractString}=nothing, cell_radius::CellRadius=nothing, spacing_factor::Real=0.95)
    @assert spacing_factor > 0 "spacing_factor must be positive. Got $(spacing_factor)."
    hex_defaults = HexagonalDefaults(cell_radius, Float64(spacing_factor))
    sink = isnothing(output) ? DataFrame : CSV.write(output)
    xml_doc = parse_file(path_to_ic_cell_xml)
    ic_cells = root(xml_doc)
    df = DataFrame(x=Float64[], y=Float64[], z=Float64[], type=String[])
    for cell_patches in child_elements(ic_cells)
        generateCellPatches!(df, cell_patches, domain_dict, hex_defaults)
    end
    free(xml_doc)

    return sink(df)
end

@deprecate generateICCell(path_to_ic_cell_xml::AbstractString, path_to_output::AbstractString, domain_dict::Dict{<:AbstractString,<:Real}) generateICCell(path_to_ic_cell_xml, domain_dict; output=path_to_output)

function generateCellPatches!(df::DataFrame, cell_patches::XMLElement, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    cell_type = attribute(cell_patches, "name")
    for patch_collection in child_elements(cell_patches)
        generatePatchCollection!(df, patch_collection, cell_type, domain_dict, hex_defaults)
    end
end

function generatePatchCollection!(df::DataFrame, patch_collection::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    patch_type = attribute(patch_collection, "type")
    for patch in child_elements(patch_collection)
        generatePatch!(patchType(patch_type), df, patch, cell_type, domain_dict, hex_defaults)
    end
end

struct DiscPatch end
struct AnnulusPatch end
struct RectanglePatch end
struct EverywherePatch end

supportedPatchTypes() = ["disc", "annulus", "rectangle", "everywhere"]

supportedPatchTextElements(::Type{DiscPatch}) = ["x0", "y0", "z0", "radius", "number", "hexagonal_grid", "spacing", "normal", "max_fails", "restrict_to_domain"]
supportedPatchTextElements(::Type{AnnulusPatch}) = ["x0", "y0", "z0", "inner_radius", "outer_radius", "number", "hexagonal_grid", "spacing", "normal", "max_fails", "restrict_to_domain"]
supportedPatchTextElements(::Type{RectanglePatch}) = ["x0", "y0", "z0", "x1", "y1", "width", "height", "number", "hexagonal_grid", "spacing", "normal", "max_fails", "restrict_to_domain"]
supportedPatchTextElements(::Type{EverywherePatch}) = ["number", "hexagonal_grid", "spacing", "normal", "max_fails"]

supportedPatchTextElements(patch_type::AbstractString) = supportedPatchTextElements(patchType(patch_type))

patchType(::Val{:disc}) = DiscPatch
patchType(::Val{:annulus}) = AnnulusPatch
patchType(::Val{:rectangle}) = RectanglePatch
patchType(::Val{:everywhere}) = EverywherePatch

function patchType(s::AbstractString)
    @assert s in supportedPatchTypes() "Patch type $(s) is not supported. Supported types are: $(supportedPatchTypes())"
    return patchType(Val(Symbol(s)))
end

function parseCenter(patch::XMLElement)
    x0 = parse(Float64, find_element(patch, "x0") |> content)
    y0 = parse(Float64, find_element(patch, "y0") |> content)
    z0 = parse(Float64, find_element(patch, "z0") |> content)
    return x0, y0, z0
end

function parseNormal(patch::XMLElement, cell_type::String)
    normal = find_element(patch, "normal")
    if isnothing(normal)
        return [0.0, 0.0, 1.0]
    else
        normal_vector = content(normal) |> x -> split(x, ",") |> x -> parse.(Float64, x)
        magnitude = sqrt(sum(normal_vector .^ 2))
        @assert magnitude > 0 "A normal vector provided for $(cell_type) with patch ID $(attribute(patch, "ID")) is 0. It must be non-zero."
        return normal_vector ./ magnitude
    end
end

function parseBoolElement(patch::XMLElement, element_name::String, default::Bool, cell_type::String)
    element = find_element(patch, element_name)
    if isnothing(element)
        return default
    end
    bool_string = strip(content(element))
    @assert bool_string ∈ ["0", "1", "true", "false"] "$(element_name) for $(cell_type) with patch ID $(attribute(patch, "ID")) must be a boolean."
    return parse(Bool, bool_string)
end

parseRestrictToDomain(patch::XMLElement, cell_type::String) = parseBoolElement(patch, "restrict_to_domain", true, cell_type)

function parseNumber(patch::XMLElement, cell_type::String)
    number = find_element(patch, "number")
    @assert !isnothing(number) "Patch ID $(attribute(patch, "ID")) for $(cell_type) must have either a number element or a hexagonal_grid element set to true."
    n = parse(Int, content(number))
    @assert n >= 0 "number for $(cell_type) with patch ID $(attribute(patch, "ID")) must be >= 0."
    return n
end

"""
    parseHexagonalGrid(patch, cell_type, hex_defaults)

Return `nothing` if the patch does not request a hexagonal grid; otherwise, return the lattice spacing (distance between neighboring cell centers).

The spacing comes from the patch's `spacing` element if present, and otherwise from `hex_defaults.spacing_factor * 2 * cell_radius`.
A patch cannot specify both `number` and `hexagonal_grid` since both determine how many cells are placed.
"""
function parseHexagonalGrid(patch::XMLElement, cell_type::String, hex_defaults::HexagonalDefaults)
    patch_descriptor = "$(cell_type) with patch ID $(attribute(patch, "ID"))"
    is_hexagonal = parseBoolElement(patch, "hexagonal_grid", false, cell_type)
    spacing_element = find_element(patch, "spacing")
    if !is_hexagonal
        @assert isnothing(spacing_element) "spacing is set for $(patch_descriptor) but hexagonal_grid is not true. spacing is only used for hexagonal grids."
        return nothing
    end
    @assert isnothing(find_element(patch, "number")) "Both number and hexagonal_grid are set for $(patch_descriptor). These are incompatible as both determine the number of cells placed in the patch."
    if !isnothing(spacing_element)
        spacing = parse(Float64, content(spacing_element))
    else
        radius = resolveCellRadius(hex_defaults.cell_radius, cell_type)
        @assert !isnothing(radius) "hexagonal_grid is true for $(patch_descriptor) but no spacing element is set and no cell_radius was passed to generateICCell for this cell type."
        spacing = hex_defaults.spacing_factor * 2 * radius
    end
    @assert spacing > 0 "spacing for $(patch_descriptor) must be positive. Got $(spacing)."
    return spacing
end

"""
    hexagonalGridPoints(spacing, half_width, half_height)

Return the `(c1, c2)` coordinates of a hexagonal lattice with the given `spacing` covering `[-half_width, half_width] x [-half_height, half_height]`.

The lattice has a point at the origin, rows parallel to the first axis separated by `spacing * sqrt(3) / 2`, and alternating rows offset by `spacing / 2`.
"""
function hexagonalGridPoints(spacing::Float64, half_width::Float64, half_height::Float64)
    row_height = spacing * sqrt(3) / 2
    n_rows = ceil(Int, half_height / row_height)
    n_cols = ceil(Int, half_width / spacing) + 1
    c1 = Float64[]
    c2 = Float64[]
    for j in -n_rows:n_rows
        y = j * row_height
        abs(y) <= half_height || continue
        offset = isodd(j) ? spacing / 2 : 0.0
        for i in -n_cols:n_cols
            x = i * spacing + offset
            abs(x) <= half_width || continue
            push!(c1, x)
            push!(c2, y)
        end
    end
    return c1, c2
end

function parseMaxFails(patch::XMLElement)
    max_fails = find_element(patch, "max_fails")
    if isnothing(max_fails)
        return 100
    end
    mf = parse(Int, content(max_fails))
    @assert mf >= 0 "Max fails must be >= 0."
    return parse(Int, content(max_fails))
end

function parseRectangleParameters(patch::XMLElement)
    x0, y0, z0 = parseCenter(patch)
    width = parseRectangleSize(patch, x0, "width", "x1")
    height = parseRectangleSize(patch, y0, "height", "y1")
    return x0, y0, z0, width, height
end

function parseRectangleSize(patch::XMLElement, c0::Float64, size_name::String, c1_name::String)
    size_element = find_element(patch, size_name)
    if !isnothing(size_element)
        return parse(Float64, content(size_element))
    end
    c1 = find_element(patch, c1_name)
    @assert !isnothing(c1) "Rectangle patch must have either a $(size_name) or $(c1_name) element."
    return parse(Float64, content(c1)) - c0
end

abstract type PatchCarveout end

struct DiscCarveout <: PatchCarveout
    x0::Float64
    y0::Float64
    z0::Float64
    radius::Float64
end

struct AnnulusCarveout <: PatchCarveout
    x0::Float64
    y0::Float64
    z0::Float64
    inner_radius::Float64
    outer_radius::Float64
end

struct RectangleCarveout <: PatchCarveout
    x0::Float64
    y0::Float64
    z0::Float64
    width::Float64
    height::Float64
end

supportedCarveoutTypes() = ["disc", "annulus", "rectangle"]

supportedCarveoutTextElements(::Type{DiscCarveout}) = ["x0", "y0", "z0", "radius"]
supportedCarveoutTextElements(::Type{AnnulusCarveout}) = ["x0", "y0", "z0", "inner_radius", "outer_radius"]
supportedCarveoutTextElements(::Type{RectangleCarveout}) = ["x0", "y0", "z0", "x1", "y1", "width", "height"]

supportedCarveoutTextElements(patch_type::AbstractString) = supportedCarveoutTextElements(carveoutType(patch_type))

carveoutType(::Val{:disc}) = DiscCarveout
carveoutType(::Val{:annulus}) = AnnulusCarveout
carveoutType(::Val{:rectangle}) = RectangleCarveout

function carveoutType(s::AbstractString)
    @assert s in supportedCarveoutTypes() "Patch type $(s) is not supported for carveouts. Supported types are: $(supportedCarveoutTypes())"
    return carveoutType(Val(Symbol(s)))
end

function parseCarveouts(patch::XMLElement)
    carveout_patches = find_element(patch, "carveout_patches")
    if isnothing(carveout_patches)
        return PatchCarveout[]
    end
    carveouts = PatchCarveout[]
    for patch_collection in child_elements(carveout_patches)
        append!(carveouts, parseCarveoutPatchCollection(patch_collection))
    end
    return carveouts
end

function parseCarveoutPatchCollection(patch_collection::XMLElement)
    patch_type = attribute(patch_collection, "type")
    @assert patch_type in ["disc", "annulus", "rectangle"] "Patch type $(patch_type) is not supported for carveouts. Supported types are: $(supportedCarveoutTypes())."
    return PatchCarveout[createCarveoutPatch(carveoutType(patch_type), patch)
                         for patch in child_elements(patch_collection)]
end

function createCarveoutPatch(::Type{DiscCarveout}, patch::XMLElement)
    x0, y0, z0 = parseCenter(patch)
    radius = parse(Float64, find_element(patch, "radius") |> content)
    return DiscCarveout(x0, y0, z0, radius)
end

function createCarveoutPatch(::Type{AnnulusCarveout}, patch::XMLElement)
    x0, y0, z0 = parseCenter(patch)
    inner_radius = parse(Float64, find_element(patch, "inner_radius") |> content)
    outer_radius = parse(Float64, find_element(patch, "outer_radius") |> content)
    return AnnulusCarveout(x0, y0, z0, inner_radius, outer_radius)
end

function createCarveoutPatch(::Type{RectangleCarveout}, patch::XMLElement)
    x0, y0, z0, width, height = parseRectangleParameters(patch)
    return RectangleCarveout(x0, y0, z0, width, height)
end

function carveOut(df::DataFrame, carveout::DiscCarveout)
    keep_ind = sqrt.((df.x .- carveout.x0) .^ 2 + (df.y .- carveout.y0) .^ 2 + (df.z .- carveout.z0) .^ 2) .> carveout.radius
    return df[keep_ind, :]
end

function carveOut(df::DataFrame, carveout::AnnulusCarveout)
    keep_ind = sqrt.((df.x .- carveout.x0) .^ 2 + (df.y .- carveout.y0) .^ 2 + (df.z .- carveout.z0) .^ 2) .> carveout.outer_radius .||
               sqrt.((df.x .- carveout.x0) .^ 2 + (df.y .- carveout.y0) .^ 2 + (df.z .- carveout.z0) .^ 2) .< carveout.inner_radius
    return df[keep_ind, :]
end

function carveOut(df::DataFrame, carveout::RectangleCarveout)
    keep_ind = (df.x .< carveout.x0) .|| (df.x .> carveout.x0 + carveout.width) .||
               (df.y .< carveout.y0) .|| (df.y .> carveout.y0 + carveout.height)
    return df[keep_ind, :]
end

function filterCells(df::DataFrame, restrict_to_domain::Bool, domain_dict::Dict{<:AbstractString,<:Real}, carveouts::Vector{PatchCarveout})
    if restrict_to_domain
        keep_ind = df.x .>= domain_dict["x_min"] .&& df.x .<= domain_dict["x_max"] .&&
                   df.y .>= domain_dict["y_min"] .&& df.y .<= domain_dict["y_max"] .&&
                   df.z .>= domain_dict["z_min"] .&& df.z .<= domain_dict["z_max"]
        df = df[keep_ind, :]
    end
    for carveout in carveouts
        df = carveOut(df, carveout)
    end
    return df
end

function createCellsDataFrame(number::Int, cell_coords_fn::Function, restrict_to_domain::Bool, domain_dict::Dict{<:AbstractString,<:Real}, carveouts::Vector{PatchCarveout}, max_fails::Int)
    total_placed = 0
    df = DataFrame(x=Float64[], y=Float64[], z=Float64[])
    fails_remaining = max_fails
    while total_placed < number
        df_new = cell_coords_fn(number - total_placed)
        df_new = filterCells(df_new, restrict_to_domain, domain_dict, carveouts)
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

"""
    annulusCoordinates(c1, c2, x0, y0, z0, normal_vector)

Map in-plane coordinates `(c1, c2)` of a disc/annulus centered at the origin to 3D coordinates of the disc/annulus centered at `(x0, y0, z0)` with the given unit `normal_vector`.
"""
function annulusCoordinates(c1::Vector{Float64}, c2::Vector{Float64}, x0::Float64, y0::Float64, z0::Float64, normal_vector::Vector{Float64})
    _df = if normal_vector[1] != 0 || normal_vector[2] != 0
        u₁ = [normal_vector[2], -normal_vector[1], 0] / sqrt(normal_vector[1]^2 + normal_vector[2]^2) # first basis vector in the plane of the disc
        u₂ = cross(normal_vector, u₁) # second basis vector in the plane of the disc
        DataFrame(x=c1 * u₁[1] + c2 * u₂[1], y=c1 * u₁[2] + c2 * u₂[2], z=c1 * u₁[3] + c2 * u₂[3])
    else
        DataFrame(x=c1, y=c2, z=zeros(length(c1)))
    end
    _df.x .+= x0
    _df.y .+= y0
    _df.z .+= z0
    return _df
end

function placeAnnulus!(df::DataFrame, inner_radius::Float64, outer_radius::Float64, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    x0, y0, z0 = parseCenter(patch)
    normal_vector = parseNormal(patch, cell_type)
    restrict_to_domain = parseRestrictToDomain(patch, cell_type)
    carveouts = parseCarveouts(patch)
    spacing = parseHexagonalGrid(patch, cell_type, hex_defaults)

    if isnothing(spacing)
        number = parseNumber(patch, cell_type)
        max_fails = parseMaxFails(patch)
        cell_coords_fn(n::Int) = begin
            r = sqrt.(inner_radius^2 .+ (outer_radius^2 - inner_radius^2) * rand(n))
            θ = 2π * rand(n)
            return annulusCoordinates(r .* cos.(θ), r .* sin.(θ), x0, y0, z0, normal_vector)
        end
        cell_df = createCellsDataFrame(number, cell_coords_fn, restrict_to_domain, domain_dict, carveouts, max_fails)
    else
        c1, c2 = hexagonalGridPoints(spacing, outer_radius, outer_radius)
        r = sqrt.(c1 .^ 2 .+ c2 .^ 2)
        keep_ind = inner_radius .<= r .<= outer_radius
        cell_df = annulusCoordinates(c1[keep_ind], c2[keep_ind], x0, y0, z0, normal_vector)
        cell_df = filterCells(cell_df, restrict_to_domain, domain_dict, carveouts)
        warnIfEmpty(cell_df, patch, cell_type)
    end
    cell_df[!, :type] .= cell_type
    append!(df, cell_df)
end

function warnIfEmpty(cell_df::DataFrame, patch::XMLElement, cell_type::String)
    if isempty(cell_df)
        @warn "Hexagonal grid for $(cell_type) with patch ID $(attribute(patch, "ID")) placed no cells. Check the spacing, patch size, domain, and carveouts."
    end
end

function generatePatch!(::Type{DiscPatch}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    radius = parse(Float64, find_element(patch, "radius") |> content)
    placeAnnulus!(df, 0.0, radius, patch, cell_type, domain_dict, hex_defaults)
end

function generatePatch!(::Type{AnnulusPatch}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    inner_radius = parse(Float64, find_element(patch, "inner_radius") |> content)
    outer_radius = parse(Float64, find_element(patch, "outer_radius") |> content)
    @assert inner_radius <= outer_radius "Inner radius must be less than or equal to outer radius. Got inner_radius=$(inner_radius) > outer_radius=$(outer_radius)."
    placeAnnulus!(df, inner_radius, outer_radius, patch, cell_type, domain_dict, hex_defaults)
end

"""
    placeBox!(df, x0, y0, z0, width, height, restrict_to_domain, patch, cell_type, domain_dict, hex_defaults)

Place cells in the axis-aligned box `[x0, x0 + width] x [y0, y0 + height]` at height `z0`, either uniformly at random or on a hexagonal grid centered in the box.
"""
function placeBox!(df::DataFrame, x0::Float64, y0::Float64, z0::Float64, width::Float64, height::Float64, restrict_to_domain::Bool, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    carveouts = parseCarveouts(patch)
    spacing = parseHexagonalGrid(patch, cell_type, hex_defaults)

    if isnothing(spacing)
        number = parseNumber(patch, cell_type)
        max_fails = parseMaxFails(patch)
        cell_coords_fn(n::Int) = DataFrame(x=x0 .+ width * rand(n), y=y0 .+ height * rand(n), z=fill(z0, n))
        cell_df = createCellsDataFrame(number, cell_coords_fn, restrict_to_domain, domain_dict, carveouts, max_fails)
    else
        c1, c2 = hexagonalGridPoints(spacing, width / 2, height / 2)
        cell_df = DataFrame(x=c1 .+ (x0 + width / 2), y=c2 .+ (y0 + height / 2), z=fill(z0, length(c1)))
        cell_df = filterCells(cell_df, restrict_to_domain, domain_dict, carveouts)
        warnIfEmpty(cell_df, patch, cell_type)
    end
    cell_df[!, :type] .= cell_type
    append!(df, cell_df)
end

function generatePatch!(::Type{RectanglePatch}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    x0, y0, z0, width, height = parseRectangleParameters(patch)
    restrict_to_domain = parseRestrictToDomain(patch, cell_type)
    placeBox!(df, x0, y0, z0, width, height, restrict_to_domain, patch, cell_type, domain_dict, hex_defaults)
end

function generatePatch!(::Type{EverywherePatch}, df::DataFrame, patch::XMLElement, cell_type::String, domain_dict::Dict{<:AbstractString,<:Real}, hex_defaults::HexagonalDefaults)
    x0 = Float64(domain_dict["x_min"])
    y0 = Float64(domain_dict["y_min"])
    width = Float64(domain_dict["x_max"] - domain_dict["x_min"])
    height = Float64(domain_dict["y_max"] - domain_dict["y_min"])
    placeBox!(df, x0, y0, 0.0, width, height, false, patch, cell_type, domain_dict, hex_defaults)
end

"""
    createICCellXMLTemplate(folder::AbstractString)

Create folder `path_to_folder` and create a template XML file for IC cells.

Introduces a new way to initialize cells in a simulation, wholly contained within a Julia package and outside a PhysiCell runtime.
It will not work in PhysiCell!
This function creates a template XML file for IC cells, showing all the current functionality of this initialization scheme.
It uses the cell type \"default\".
Create ICs for more cell types by copying the `cell_patches` element.
"""
function createICCellXMLTemplate(path_to_folder::AbstractString)
    path_to_ic_cell_xml = joinpath(path_to_folder, "cells.xml")
    mkpath(dirname(path_to_ic_cell_xml))
    xml_doc = XMLDocument()
    xml_root = create_root(xml_doc, "ic_cells")

    e_patches = new_child(xml_root, "cell_patches")
    set_attribute(e_patches, "name", "default")

    # make disc patch
    e_discs = new_child(e_patches, "patch_collection")
    set_attribute(e_discs, "type", "disc")
    e_patch = new_child(e_discs, "patch")
    set_attribute(e_patch, "ID", "1")
    for (name, value) in [("x0", "0.0"), ("y0", "0.0"), ("z0", "0.0"), ("radius", "40.0"), ("number", "50"), ("normal", "0,0,1"), ("max_fails", "100")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    ## make rectangle carveouts for the disc
    e_carveout_patches = new_child(e_patch, "carveout_patches")
    e_rectangle_carveouts = new_child(e_carveout_patches, "patch_collection")
    set_attribute(e_rectangle_carveouts, "type", "rectangle")
    e_carveout_patch = new_child(e_rectangle_carveouts, "patch")
    set_attribute(e_carveout_patch, "ID", "1")
    for (name, value) in [("x0", "0.0"), ("y0", "0.0"), ("z0", "0.0"), ("width", "10.0"), ("height", "10.0")]
        e = new_child(e_carveout_patch, name)
        set_content(e, value)
    end

    ## make a second disc patch filled with a hexagonal grid of cells
    e_patch = new_child(e_discs, "patch")
    set_attribute(e_patch, "ID", "2")
    for (name, value) in [("x0", "200.0"), ("y0", "200.0"), ("z0", "0.0"), ("radius", "40.0"), ("hexagonal_grid", "true"), ("spacing", "12.0")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    # make annulus patch
    e_annuli = new_child(e_patches, "patch_collection")
    set_attribute(e_annuli, "type", "annulus")
    e_patch = new_child(e_annuli, "patch")
    set_attribute(e_patch, "ID", "1")
    for (name, value) in [("x0", "50.0"), ("y0", "50.0"), ("z0", "0.0"), ("inner_radius", "10.0"), ("outer_radius", "200.0"), ("number", "50"), ("restrict_to_domain", "true")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    # make rectangle patch
    e_rectangles = new_child(e_patches, "patch_collection")
    set_attribute(e_rectangles, "type", "rectangle")
    e_patch = new_child(e_rectangles, "patch")
    set_attribute(e_patch, "ID", "1")
    for (name, value) in [("x0", "-50.0"), ("y0", "-50.0"), ("z0", "0.0"), ("x1", "50.0"), ("y1", "50.0"), ("number", "10")]
        e = new_child(e_patch, name)
        set_content(e, value)
    end

    ## make disc carveouts for the rectangle
    e_carveout_patches = new_child(e_patch, "carveout_patches")
    e_disc_carveouts = new_child(e_carveout_patches, "patch_collection")
    set_attribute(e_disc_carveouts, "type", "disc")
    e_carveout_patch = new_child(e_disc_carveouts, "patch")
    set_attribute(e_carveout_patch, "ID", "1")
    for (name, value) in [("x0", "0.0"), ("y0", "0.0"), ("z0", "0.0"), ("radius", "10.0")]
        e = new_child(e_carveout_patch, name)
        set_content(e, value)
    end

    # make everywhere patch
    e_everywhere = new_child(e_patches, "patch_collection")
    set_attribute(e_everywhere, "type", "everywhere")
    e_patch = new_child(e_everywhere, "patch")
    set_attribute(e_patch, "ID", "1")
    e = new_child(e_patch, "number")
    set_content(e, "78")

    ## make annulus carveouts for the everywhere patch
    e_carveout_patches = new_child(e_patch, "carveout_patches")
    e_annuli_carveouts = new_child(e_carveout_patches, "patch_collection")
    set_attribute(e_annuli_carveouts, "type", "annulus")
    e_carveout_patch = new_child(e_annuli_carveouts, "patch")
    set_attribute(e_carveout_patch, "ID", "1")
    for (name, value) in [("x0", "0.0"), ("y0", "0.0"), ("z0", "0.0"), ("inner_radius", "100.0"), ("outer_radius", "400.0")]
        e = new_child(e_carveout_patch, name)
        set_content(e, value)
    end

    save_file(xml_doc, path_to_ic_cell_xml)
    free(xml_doc)
end

end
