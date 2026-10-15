module PhysiCellCellCreator

using LightXML, CSV, DataFrames, LinearAlgebra, Compat

export generateICCell, createICCellXMLTemplate

@compat public supportedPatchTypes, supportedPatchTextElements, supportedPlacements, supportedCarveoutTypes, supportedCarveoutTextElements

include("frame.jl")
include("shapes.jl")
include("parse.jl")
include("carveouts.jl")
include("patches.jl")
include("template.jl")

"""
    generateICCell(path_to_ic_cell_xml, domain_dict; output=nothing, cell_radius=nothing, spacing_factor=0.95, boundary_tolerance=1e-10)

Generate initial cell positions from the IC cell XML file at `path_to_ic_cell_xml`.

`domain_dict` must contain the keys `"x_min"`, `"x_max"`, `"y_min"`, `"y_max"`, `"z_min"`, and `"z_max"`.

# Keyword arguments
- `output`: if `nothing` (default), return a `DataFrame` with columns `x`, `y`, `z`, and `type`; otherwise, write a CSV file to this path.
- `cell_radius`: the radius of the cells being placed, used to compute the default spacing for patches that use a hexagonal grid but do not specify a `spacing` element.
  Either a single `Real` applied to every cell type or a `Dict` mapping cell type names (the `name` attribute of `cell_patches`) to radii.
  Only required when a hexagonal grid patch omits `spacing`.
- `spacing_factor`: multiplies the cell diameter to get the default spacing, i.e., `spacing = spacing_factor * 2 * cell_radius`.
  The default of `0.95` matches the slightly compressed packing used in PhysiCell's sample projects; `1.0` places cells tangent to their neighbors.
- `boundary_tolerance`: for hexagonal grids, lattice points within `boundary_tolerance * spacing` outside the patch boundary are still kept.
  The default of `1e-10` only absorbs floating-point rounding so that points exactly on the boundary are not dropped; set it to `0` for a strict comparison.
"""
function generateICCell(path_to_ic_cell_xml::AbstractString, domain_dict::Dict{<:AbstractString,<:Real}; output::Union{Nothing,AbstractString}=nothing, cell_radius::CellRadius=nothing, spacing_factor::Real=0.95, boundary_tolerance::Real=1e-10)
    @assert spacing_factor > 0 "spacing_factor must be positive. Got $(spacing_factor)."
    @assert boundary_tolerance >= 0 "boundary_tolerance must be non-negative. Got $(boundary_tolerance)."
    hex_defaults = HexagonalDefaults(cell_radius, Float64(spacing_factor), Float64(boundary_tolerance))
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

end
