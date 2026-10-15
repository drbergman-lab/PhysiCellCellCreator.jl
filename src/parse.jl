patchDescriptor(patch::XMLElement, cell_type::String) = "$(cell_type) with patch ID $(attribute(patch, "ID"))"

function findRequiredElement(patch::XMLElement, name::String, descriptor::String)
    element = find_element(patch, name)
    @assert !isnothing(element) "$(descriptor) must have a $(name) element."
    return element
end

parseFloatElement(patch::XMLElement, name::String, descriptor::String) = parse(Float64, strip(content(findRequiredElement(patch, name, descriptor))))

function parseOptionalFloatElement(patch::XMLElement, name::String, default::Float64)
    element = find_element(patch, name)
    return isnothing(element) ? default : parse(Float64, strip(content(element)))
end

function parsePositiveFloatElement(patch::XMLElement, name::String, descriptor::String)
    value = parseFloatElement(patch, name, descriptor)
    @assert value > 0 "$(name) for $(descriptor) must be positive. Got $(value)."
    return value
end

parseCenter(patch::XMLElement, descriptor::String) = [parseFloatElement(patch, name, descriptor) for name in ["x0", "y0", "z0"]]

function parseNormal(patch::XMLElement, descriptor::String)
    normal = find_element(patch, "normal")
    if isnothing(normal)
        return [0.0, 0.0, 1.0]
    end
    normal_vector = content(normal) |> x -> split(x, ",") |> x -> parse.(Float64, strip.(x))
    @assert length(normal_vector) == 3 "The normal vector for $(descriptor) must have three comma-separated components. Got $(content(normal))."
    magnitude = sqrt(sum(normal_vector .^ 2))
    @assert magnitude > 0 "A normal vector provided for $(descriptor) is 0. It must be non-zero."
    return normal_vector ./ magnitude
end

const pi_multiple_regex = r"^([+-]?)\s*(\d+\.?\d*|\.\d+)?\s*\*?\s*(?:pi|π)\s*(?:/\s*(\d+\.?\d*|\.\d+))?$"

"""
    parsePiMultiple(s, descriptor)

Parse a rotation written as a multiple of pi, e.g. `pi`, `-pi/2`, `2pi/3`, `2*pi/3`, `0.5pi`, or `π/6`, into radians without evaluating it as code.
"""
function parsePiMultiple(s::AbstractString, descriptor::String)
    m = match(pi_multiple_regex, s)
    @assert !isnothing(m) """
    The rotation for $(descriptor) contains pi but could not be parsed: $(s)
    Accepted forms are a (signed) coefficient, pi or π, and an optional divisor, e.g.
        <rotation>pi/6</rotation>, <rotation>-2pi/3</rotation>, <rotation>0.5*pi</rotation>
    """
    sign = m[1] == "-" ? -1.0 : 1.0
    coefficient = isnothing(m[2]) ? 1.0 : parse(Float64, m[2])
    denominator = isnothing(m[3]) ? 1.0 : parse(Float64, m[3])
    @assert denominator != 0 "The rotation for $(descriptor) divides by zero: $(s)"
    return sign * coefficient * π / denominator
end

"""
    parseRotation(patch, descriptor)

Parse the optional `rotation` element into radians, following the rules of PhysiCellECMCreator.jl.

The value is in degrees by default and in radians if it contains `pi` or `π`.
A `units` attribute starting with `rad` or `deg` sets the units explicitly; using `pi` with degrees is an error.
"""
function parseRotation(patch::XMLElement, descriptor::String)
    rotation_element = find_element(patch, "rotation")
    if isnothing(rotation_element)
        return 0.0
    end
    θ = strip(content(rotation_element))
    units = attribute(rotation_element, "units")
    uses_pi = occursin("pi", θ) || occursin("π", θ)
    if isnothing(units)
        is_radians = uses_pi
    else
        units_lower = lowercase(strip(units))
        @assert startswith(units_lower, "rad") || startswith(units_lower, "deg") "The units attribute of the rotation for $(descriptor) must start with \"rad\" or \"deg\". Got \"$(units)\"."
        is_radians = startswith(units_lower, "rad")
    end
    @assert !uses_pi || is_radians """
    The rotation for $(descriptor) uses pi or π, so its units cannot be degrees. Solutions:
    - Omit the units attribute: <rotation>pi/6</rotation> # interpreted as radians
    - Specify the units as radians: <rotation units="rad">pi/6</rotation>
    - Convert to degrees: <rotation units="deg">30</rotation> OR <rotation>30</rotation> # interpreted as degrees since it does not contain pi or π
    """
    value = uses_pi ? parsePiMultiple(θ, descriptor) : parse(Float64, θ)
    return is_radians ? value : value * π / 180
end

function parseBoolElement(patch::XMLElement, element_name::String, default::Bool, descriptor::String)
    element = find_element(patch, element_name)
    if isnothing(element)
        return default
    end
    bool_string = strip(content(element))
    @assert bool_string ∈ ["0", "1", "true", "false"] "$(element_name) for $(descriptor) must be a boolean."
    return parse(Bool, bool_string)
end

parseRestrictToDomain(patch::XMLElement, descriptor::String) = parseBoolElement(patch, "restrict_to_domain", true, descriptor)

function parseNumber(patch::XMLElement, descriptor::String)
    number = find_element(patch, "number")
    @assert !isnothing(number) "$(descriptor) uses random placement and so must have a number element."
    n = parse(Int, strip(content(number)))
    @assert n >= 0 "number for $(descriptor) must be >= 0."
    return n
end

function parseMaxFails(patch::XMLElement, descriptor::String)
    max_fails = find_element(patch, "max_fails")
    if isnothing(max_fails)
        return 100
    end
    mf = parse(Int, strip(content(max_fails)))
    @assert mf >= 0 "max_fails for $(descriptor) must be >= 0."
    return mf
end

"""
    parseRectangleParameters(patch, descriptor)

Return the center, width, and height of a rectangle given by its minimum corner `(x0, y0)` at height `z0` and either `width`/`height` or the opposite corner `x1`/`y1`.
"""
function parseRectangleParameters(patch::XMLElement, descriptor::String)
    x0, y0, z0 = parseCenter(patch, descriptor)
    width = parseRectangleSize(patch, x0, "width", "x1", descriptor)
    height = parseRectangleSize(patch, y0, "height", "y1", descriptor)
    return [x0 + width / 2, y0 + height / 2, z0], width, height
end

function parseRectangleSize(patch::XMLElement, c0::Float64, size_name::String, c1_name::String, descriptor::String)
    size_element = find_element(patch, size_name)
    size = if !isnothing(size_element)
        parse(Float64, strip(content(size_element)))
    else
        c1 = find_element(patch, c1_name)
        @assert !isnothing(c1) "Rectangle $(descriptor) must have either a $(size_name) or $(c1_name) element."
        parse(Float64, strip(content(c1))) - c0
    end
    @assert size > 0 "Rectangle $(descriptor) must have positive $(size_name). Got $(size) (if using $(c1_name), it must be greater than the corresponding corner coordinate)."
    return size
end
