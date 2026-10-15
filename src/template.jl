function addPatch!(patch_collection::XMLElement, id::String, elements::Vector{Pair{String,String}}; placement::Union{Nothing,String}=nothing)
    e_patch = new_child(patch_collection, "patch")
    set_attribute(e_patch, "ID", id)
    if !isnothing(placement)
        set_attribute(e_patch, "placement", placement)
    end
    for (name, value) in elements
        e = new_child(e_patch, name)
        set_content(e, value)
    end
    return e_patch
end

function addPatchCollection!(parent::XMLElement, patch_type::String)
    e_collection = new_child(parent, "patch_collection")
    set_attribute(e_collection, "type", patch_type)
    return e_collection
end

function addCarveout!(e_patch::XMLElement, carveout_type::String, id::String, elements::Vector{Pair{String,String}})
    e_carveout_patches = find_element(e_patch, "carveout_patches")
    if isnothing(e_carveout_patches)
        e_carveout_patches = new_child(e_patch, "carveout_patches")
    end
    e_collection = addPatchCollection!(e_carveout_patches, carveout_type)
    return addPatch!(e_collection, id, elements)
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

    # disc patches: one random with a rectangle carved out, one filled with a hexagonal grid
    e_discs = addPatchCollection!(e_patches, "disc")
    e_patch = addPatch!(e_discs, "1", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "50", "normal" => "0,0,1", "max_fails" => "100"])
    addCarveout!(e_patch, "rectangle", "1", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "width" => "10.0", "height" => "10.0"])
    addPatch!(e_discs, "2", ["x0" => "200.0", "y0" => "200.0", "z0" => "0.0", "radius" => "40.0", "spacing" => "12.0"]; placement="hexagonal")

    # annulus patch
    e_annuli = addPatchCollection!(e_patches, "annulus")
    addPatch!(e_annuli, "1", ["x0" => "50.0", "y0" => "50.0", "z0" => "0.0", "inner_radius" => "10.0", "outer_radius" => "200.0", "number" => "50", "restrict_to_domain" => "true"])

    # rectangle patch with a disc carved out
    e_rectangles = addPatchCollection!(e_patches, "rectangle")
    e_patch = addPatch!(e_rectangles, "1", ["x0" => "-50.0", "y0" => "-50.0", "z0" => "0.0", "x1" => "50.0", "y1" => "50.0", "number" => "10"])
    addCarveout!(e_patch, "disc", "1", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "10.0"])

    # elliptical disc patch (random placement) with a shell of hexagonally placed cells around it (ellipse patch)
    e_elliptical_discs = addPatchCollection!(e_patches, "elliptical_disc")
    addPatch!(e_elliptical_discs, "1", ["x0" => "-200.0", "y0" => "200.0", "z0" => "0.0", "a" => "80.0", "b" => "40.0", "rotation" => "pi/6", "number" => "60"])
    e_ellipses = addPatchCollection!(e_patches, "ellipse")
    addPatch!(e_ellipses, "1", ["x0" => "-200.0", "y0" => "200.0", "z0" => "0.0", "a" => "80.0", "b" => "40.0", "thickness" => "20.0", "rotation" => "pi/6", "spacing" => "10.0"]; placement="hexagonal")

    # everywhere patch avoiding an annulus and the ellipse (with its shell) above
    e_everywhere = addPatchCollection!(e_patches, "everywhere")
    e_patch = addPatch!(e_everywhere, "1", ["number" => "78"])
    addCarveout!(e_patch, "annulus", "1", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "100.0", "outer_radius" => "400.0"])
    addCarveout!(e_patch, "elliptical_disc", "1", ["x0" => "-200.0", "y0" => "200.0", "z0" => "0.0", "a" => "100.0", "b" => "60.0", "rotation" => "pi/6"])

    save_file(xml_doc, path_to_ic_cell_xml)
    free(xml_doc)
end
