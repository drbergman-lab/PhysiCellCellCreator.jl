"""
    Carveout(shape, frame)

A region from which cells are removed after placement.
A planar carveout removes every cell whose projection onto the carveout's plane lies inside the shape, i.e. it is an infinite prism along the carveout's normal.
"""
struct Carveout{S<:Shape}
    shape::S
    frame::Frame
end

"""
    supportedCarveoutTypes()

The values accepted by the `type` attribute of a `patch_collection` inside a `carveout_patches` element: every patch type except `everywhere`, since an `everywhere` carveout would remove every cell of its patch.
"""
supportedCarveoutTypes() = filter(!=("everywhere"), supportedPatchTypes())

function carveoutType(s::AbstractString)
    @assert s in supportedCarveoutTypes() "Patch type $(s) is not supported for carveouts. Supported types are: $(supportedCarveoutTypes())"
    return patchType(s)
end

"""
    supportedCarveoutTextElements(carveout_type)

The elements a carveout `patch` of the given type may contain: the geometric elements of the patch type of the same name.
"""
supportedCarveoutTextElements(carveout_type::AbstractString) = geometryElements(carveoutType(carveout_type))

function parseCarveouts(patch::XMLElement, descriptor::String)
    carveout_patches = find_element(patch, "carveout_patches")
    if isnothing(carveout_patches)
        return Carveout[]
    end
    carveouts = Carveout[]
    for patch_collection in child_elements(carveout_patches)
        T = carveoutType(attribute(patch_collection, "type"))
        for carveout_patch in child_elements(patch_collection)
            carveout_descriptor = "carveout ID $(attribute(carveout_patch, "ID")) of $(descriptor)"
            shape, frame = parseGeometry(T, carveout_patch, carveout_descriptor)
            push!(carveouts, Carveout(shape, frame))
        end
    end
    return carveouts
end

function carveOut!(df::DataFrame, carveout::Carveout{<:PlanarShape})
    c1, c2, _ = toLocal(carveout.frame, df)
    deleteat!(df, findall(inside(carveout.shape, c1, c2)))
end

"""
    filterCells!(df, restrict_to_domain, domain_dict, carveouts)

Remove in place the rows of `df` that fall outside the domain (if `restrict_to_domain`) or inside any of the `carveouts`.
"""
function filterCells!(df::DataFrame, restrict_to_domain::Bool, domain_dict::Dict{<:AbstractString,<:Real}, carveouts::Vector{Carveout})
    if restrict_to_domain
        keep_ind = df.x .>= domain_dict["x_min"] .&& df.x .<= domain_dict["x_max"] .&&
                   df.y .>= domain_dict["y_min"] .&& df.y .<= domain_dict["y_max"] .&&
                   df.z .>= domain_dict["z_min"] .&& df.z .<= domain_dict["z_max"]
        deleteat!(df, findall(.!keep_ind))
    end
    for carveout in carveouts
        carveOut!(df, carveout)
    end
    return df
end
