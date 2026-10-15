using LightXML, DataFrames, LinearAlgebra, Random

Random.seed!(20261014)

shape_folder = "shape_test_folder"
mkpath(shape_folder)
shape_xml = joinpath(shape_folder, "cells.xml")
shape_domain = Dict{String,Float64}("x_min" => -500.0, "x_max" => 500.0, "y_min" => -500.0, "y_max" => 500.0, "z_min" => -500.0, "z_max" => 500.0)

function writeShapeXML(path::String, patch_type::String, elements::Vector{Pair{String,String}}; placement::Union{Nothing,String}=nothing, carveouts::Vector=[])
    xml_doc = XMLDocument()
    xml_root = create_root(xml_doc, "ic_cells")
    e_patches = new_child(xml_root, "cell_patches")
    set_attribute(e_patches, "name", "default")
    e_collection = new_child(e_patches, "patch_collection")
    set_attribute(e_collection, "type", patch_type)
    e_patch = new_child(e_collection, "patch")
    set_attribute(e_patch, "ID", "1")
    isnothing(placement) || set_attribute(e_patch, "placement", placement)
    for (name, value) in elements
        set_content(new_child(e_patch, name), value)
    end
    if !isempty(carveouts)
        e_carveouts = new_child(e_patch, "carveout_patches")
        for (i, (carveout_type, carveout_elements)) in enumerate(carveouts)
            e_carveout_collection = new_child(e_carveouts, "patch_collection")
            set_attribute(e_carveout_collection, "type", carveout_type)
            e_carveout = new_child(e_carveout_collection, "patch")
            set_attribute(e_carveout, "ID", string(i))
            for (name, value) in carveout_elements
                set_content(new_child(e_carveout, name), value)
            end
        end
    end
    save_file(xml_doc, path)
    free(xml_doc)
    return path
end

# angle (mod π) of the principal axis of a 2D point cloud
function principalAxisAngle(x::AbstractVector, y::AbstractVector)
    dx = x .- sum(x) / length(x)
    dy = y .- sum(y) / length(y)
    sxx, syy, sxy = sum(dx .^ 2), sum(dy .^ 2), sum(dx .* dy)
    return 0.5 * atan(2 * sxy, sxx - syy)
end

angleDifference(θ1, θ2) = abs(rem(θ1 - θ2, π, RoundNearest))

ellipseFunctional(c1, c2, a, b) = (c1 / a)^2 + (c2 / b)^2

# local in-plane coordinates of (x, y) relative to a center rotated by θ counter-clockwise
function rotateBack(x, y, x0, y0, θ)
    dx, dy = x .- x0, y .- y0
    return cos(θ) .* dx .+ sin(θ) .* dy, -sin(θ) .* dx .+ cos(θ) .* dy
end

function nearestNeighborDistance(df::DataFrame)
    n = size(df, 1)
    dists = fill(Inf, n)
    for i in 1:n, j in 1:n
        i == j && continue
        dists[i] = min(dists[i], sqrt((df.x[i] - df.x[j])^2 + (df.y[i] - df.y[j])^2 + (df.z[i] - df.z[j])^2))
    end
    return dists
end

@testset "elliptical disc: uniform and inside" begin
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "10.0", "y0" => "-20.0", "z0" => "3.0", "a" => "120.0", "b" => "60.0", "number" => "100000"])
    df = generateICCell(shape_xml, shape_domain)
    @test size(df, 1) == 100000
    @test all(df.z .== 3.0)
    f = ellipseFunctional.(df.x .- 10.0, df.y .+ 20.0, 120.0, 60.0)
    @test all(f .<= 1 + 1e-12)
    # uniform in area: the ellipse scaled by s holds a fraction s^2 of the cells
    @test isapprox(count(f .<= 0.5^2) / size(df, 1), 0.25; atol=0.01)
    @test isapprox(count(f .<= 0.8^2) / size(df, 1), 0.64; atol=0.01)
    @test isapprox(count(df.x .> 10.0) / size(df, 1), 0.5; atol=0.01)
    @test angleDifference(principalAxisAngle(df.x, df.y), 0.0) < deg2rad(1)
end

@testset "elliptical disc: rotation" begin
    for (rotation, units) in [("pi/6", nothing), ("30", nothing), ("30", "deg"), ("π/6", "rad"), ("0.5235987755982988", "radians")]
        elements = ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "number" => "50000"]
        writeShapeXML(shape_xml, "elliptical_disc", elements)
        # add the rotation element by hand so the units attribute can be set
        xml_doc = parse_file(shape_xml)
        e_patch = find_element(find_element(find_element(root(xml_doc), "cell_patches"), "patch_collection"), "patch")
        e_rotation = new_child(e_patch, "rotation")
        set_content(e_rotation, rotation)
        isnothing(units) || set_attribute(e_rotation, "units", units)
        save_file(xml_doc, shape_xml)
        free(xml_doc)

        df = generateICCell(shape_xml, shape_domain)
        c1, c2 = rotateBack(df.x, df.y, 0.0, 0.0, π / 6)
        @test all(ellipseFunctional.(c1, c2, 120.0, 40.0) .<= 1 + 1e-12)
        @test angleDifference(principalAxisAngle(df.x, df.y), π / 6) < deg2rad(1)
    end
    # a negative rotation turns the other way
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "rotation" => "-pi/3", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test angleDifference(principalAxisAngle(df.x, df.y), -π / 3) < deg2rad(1)
end

@testset "elliptical disc: normal and frame continuity" begin
    # normal along y: the disc lies in the plane y = y0 with a along x
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "25.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "normal" => "0,1,0", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(abs.(df.y .- 25.0) .< 1e-10)
    @test all(ellipseFunctional.(df.x, df.z, 120.0, 40.0) .<= 1 + 1e-12)
    @test angleDifference(principalAxisAngle(df.x, df.z), 0.0) < deg2rad(1)

    # a barely tilted normal gives a barely tilted frame: a still runs along x
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "normal" => "1e-6,0,1", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test angleDifference(principalAxisAngle(df.x, df.y), 0.0) < deg2rad(1)
    n = [1e-6, 0.0, 1.0] / norm([1e-6, 0.0, 1.0])
    @test all(abs.(df.x .* n[1] .+ df.y .* n[2] .+ df.z .* n[3]) .< 1e-9) # in the plane

    # a normal nearly parallel to x still uses the projected x-axis (almost -z) for a, not the y-axis fallback
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "normal" => "1,0,1e-6", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(abs.(df.x .+ 1e-6 .* df.z) .< 1e-9) # in the plane
    @test angleDifference(principalAxisAngle(df.y, df.z), π / 2) < deg2rad(1)

    # normal along x: a runs along y (the projection of the y-axis), rotation turns about x
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "normal" => "2,0,0", "rotation" => "pi/4", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(abs.(df.x) .< 1e-10)
    @test angleDifference(principalAxisAngle(df.y, df.z), π / 4) < deg2rad(1)

    # malformed normals
    for normal in ["0,0,0", "1,0"]
        writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "120.0", "b" => "40.0", "normal" => normal, "number" => "10"])
        @test_throws AssertionError generateICCell(shape_xml, shape_domain)
    end
end

@testset "ellipse: ring with thickness" begin
    writeShapeXML(shape_xml, "ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "thickness" => "30.0", "rotation" => "pi/3", "number" => "50000"])
    df = generateICCell(shape_xml, shape_domain)
    @test size(df, 1) == 50000
    c1, c2 = rotateBack(df.x, df.y, 0.0, 0.0, π / 3)
    @test all(ellipseFunctional.(c1, c2, 130.0, 80.0) .<= 1 + 1e-12)
    @test all(ellipseFunctional.(c1, c2, 100.0, 50.0) .>= 1 - 1e-12)
    @test angleDifference(principalAxisAngle(df.x, df.y), π / 3) < deg2rad(1)
    # uniform in area: the fraction inside the ellipse with axes (115, 65) is the area ratio
    expected = (115 * 65 - 100 * 50) / (130 * 80 - 100 * 50)
    @test isapprox(count(ellipseFunctional.(c1, c2, 115.0, 65.0) .<= 1) / size(df, 1), expected; atol=0.01)

    # a very thin ring is sampled in bounded batches
    writeShapeXML(shape_xml, "ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "thickness" => "0.01", "number" => "200"])
    df = generateICCell(shape_xml, shape_domain)
    @test size(df, 1) == 200
    @test all(ellipseFunctional.(df.x, df.y, 100.01, 50.01) .<= 1 + 1e-12)
    @test all(ellipseFunctional.(df.x, df.y, 100.0, 50.0) .>= 1 - 1e-12)
    @test PhysiCellCellCreator.max_proposal_batch == 1_000_000

    # a tolerance larger than the inner axes keeps everything inside the (enlarged) outer ellipse
    ring = PhysiCellCellCreator.EllipticalRing(100.0, 50.0, 30.0)
    @test PhysiCellCellCreator.inside(ring, [0.0, 140.0, 200.0], [0.0, 0.0, 0.0]; tol=60.0) == [true, true, false]
    @test PhysiCellCellCreator.inside(ring, [0.0, 140.0], [0.0, 0.0]; tol=0.0) == [false, false]

    # thickness is required and must be positive
    writeShapeXML(shape_xml, "ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "number" => "10"])
    @test_throws AssertionError generateICCell(shape_xml, shape_domain)
    writeShapeXML(shape_xml, "ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "thickness" => "0.0", "number" => "10"])
    @test_throws AssertionError generateICCell(shape_xml, shape_domain)
end

@testset "disc and annulus regressions" begin
    writeShapeXML(shape_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "5.0", "radius" => "40.0", "number" => "40000"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(df.z .== 5.0)
    r = sqrt.(df.x .^ 2 .+ df.y .^ 2)
    @test all(r .<= 40.0 + 1e-12)
    @test isapprox(count(r .<= 20.0) / size(df, 1), 0.25; atol=0.01)

    writeShapeXML(shape_xml, "annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "20.0", "outer_radius" => "60.0", "number" => "40000"])
    df = generateICCell(shape_xml, shape_domain)
    r = sqrt.(df.x .^ 2 .+ df.y .^ 2)
    @test all(20.0 - 1e-12 .<= r .<= 60.0 + 1e-12)
    @test isapprox(count(r .<= 40.0) / size(df, 1), (40^2 - 20^2) / (60^2 - 20^2); atol=0.01)

    # inner_radius of 0 is a disc; inner_radius equal to outer_radius is a circle
    writeShapeXML(shape_xml, "annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "0.0", "outer_radius" => "60.0", "number" => "100"])
    @test size(generateICCell(shape_xml, shape_domain), 1) == 100
    writeShapeXML(shape_xml, "annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "60.0", "outer_radius" => "60.0", "number" => "100"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(isapprox.(sqrt.(df.x .^ 2 .+ df.y .^ 2), 60.0; atol=1e-9))
    writeShapeXML(shape_xml, "annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "70.0", "outer_radius" => "60.0", "number" => "100"])
    @test_throws AssertionError generateICCell(shape_xml, shape_domain)

    # a rotation on a disc only turns its hexagonal lattice
    writeShapeXML(shape_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "30.0", "spacing" => "10.0", "rotation" => "90"]; placement="hexagonal")
    df = generateICCell(shape_xml, shape_domain)
    hasCellAt(x, y) = any(isapprox.(df.x, x; atol=1e-8) .&& isapprox.(df.y, y; atol=1e-8))
    @test hasCellAt(0.0, 10.0)
    @test !hasCellAt(10.0, 0.0)
end

@testset "rectangle regressions" begin
    writeShapeXML(shape_xml, "rectangle", ["x0" => "100.0", "y0" => "-50.0", "z0" => "7.0", "width" => "50.0", "height" => "100.0", "number" => "40000"])
    df = generateICCell(shape_xml, shape_domain)
    @test all(100.0 .<= df.x .<= 150.0)
    @test all(-50.0 .<= df.y .<= 50.0)
    @test all(df.z .== 7.0)
    @test isapprox(count(df.x .< 125.0) / size(df, 1), 0.5; atol=0.01)
    @test isapprox(count(df.y .< -25.0) / size(df, 1), 0.25; atol=0.01)
end

@testset "elliptical carveouts" begin
    # a rotated elliptical disc carved out of a rectangle: no cell inside it, and the right number removed from a lattice
    carveout = ("elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "rotation" => "pi/4"])
    writeShapeXML(shape_xml, "rectangle", ["x0" => "-200.0", "y0" => "-200.0", "z0" => "0.0", "width" => "400.0", "height" => "400.0", "number" => "20000"]; carveouts=[carveout])
    df = generateICCell(shape_xml, shape_domain)
    @test size(df, 1) == 20000
    c1, c2 = rotateBack(df.x, df.y, 0.0, 0.0, π / 4)
    @test all(ellipseFunctional.(c1, c2, 100.0, 50.0) .> 1)

    spacing = 4.0
    writeShapeXML(shape_xml, "rectangle", ["x0" => "-200.0", "y0" => "-200.0", "z0" => "0.0", "width" => "400.0", "height" => "400.0", "spacing" => string(spacing)]; placement="hexagonal")
    n_full = size(generateICCell(shape_xml, shape_domain), 1)
    writeShapeXML(shape_xml, "rectangle", ["x0" => "-200.0", "y0" => "-200.0", "z0" => "0.0", "width" => "400.0", "height" => "400.0", "spacing" => string(spacing)]; placement="hexagonal", carveouts=[carveout])
    n_carved = size(generateICCell(shape_xml, shape_domain), 1)
    expected_removed = π * 100 * 50 / (spacing^2 * sqrt(3) / 2)
    @test isapprox(n_full - n_carved, expected_removed; rtol=0.03)

    # an ellipse (ring) carveout leaves the inside of the inner ellipse and the outside of the outer one
    carveout = ("ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "thickness" => "30.0"])
    writeShapeXML(shape_xml, "everywhere", ["number" => "20000"]; carveouts=[carveout])
    df = generateICCell(shape_xml, shape_domain)
    f_inner = ellipseFunctional.(df.x, df.y, 100.0, 50.0)
    f_outer = ellipseFunctional.(df.x, df.y, 130.0, 80.0)
    @test all(f_inner .< 1 .|| f_outer .> 1)
    @test any(f_inner .< 1)
    @test any(f_outer .> 1)
end

@testset "carveouts are prisms along their normal" begin
    # cells in the plane x = 0; a disc carveout with the default normal removes the band |y| <= 20 for every z
    carveout = ("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "20.0"])
    writeShapeXML(shape_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "60.0", "normal" => "1,0,0", "number" => "5000"]; carveouts=[carveout])
    df = generateICCell(shape_xml, shape_domain)
    @test all(abs.(df.y) .> 20.0)
    @test any(abs.(df.y) .> 20.0 .&& abs.(df.z) .< 5.0)
    @test any(sqrt.(df.y .^ 2 .+ df.z .^ 2) .< 30.0) # cells near the carveout that a ball would have removed survive

    # a carveout with its own normal: a disc normal to x removes cells with sqrt(y^2 + z^2) <= 20 from a patch in the xy-plane
    carveout = ("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "20.0", "normal" => "1,0,0"])
    writeShapeXML(shape_xml, "rectangle", ["x0" => "-100.0", "y0" => "-100.0", "z0" => "0.0", "width" => "200.0", "height" => "200.0", "number" => "5000"]; carveouts=[carveout])
    df = generateICCell(shape_xml, shape_domain)
    @test all(abs.(df.y) .> 20.0)
    @test any(abs.(df.x) .< 20.0)

    # a carveout at a different height still applies (prism), where a ball would not
    carveout = ("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "300.0", "radius" => "20.0"])
    writeShapeXML(shape_xml, "rectangle", ["x0" => "-100.0", "y0" => "-100.0", "z0" => "0.0", "width" => "200.0", "height" => "200.0", "number" => "5000"]; carveouts=[carveout])
    df = generateICCell(shape_xml, shape_domain)
    @test all(sqrt.(df.x .^ 2 .+ df.y .^ 2) .> 20.0)
end

@testset "hexagonal elliptical patches" begin
    spacing = 4.0
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "spacing" => string(spacing)]; placement="hexagonal")
    df = generateICCell(shape_xml, shape_domain)
    @test all(ellipseFunctional.(df.x, df.y, 100.0, 50.0) .<= 1 + 1e-9)
    @test any(abs.(df.x) .< 1e-8 .&& abs.(df.y) .< 1e-8)
    @test all(isapprox.(nearestNeighborDistance(df), spacing; atol=1e-8))
    @test isapprox(size(df, 1), π * 100 * 50 / (spacing^2 * sqrt(3) / 2); rtol=0.03)

    # rotation turns the ellipse and its lattice together
    writeShapeXML(shape_xml, "elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "rotation" => "pi/6", "spacing" => "10.0"]; placement="hexagonal")
    df = generateICCell(shape_xml, shape_domain)
    c1, c2 = rotateBack(df.x, df.y, 0.0, 0.0, π / 6)
    @test all(ellipseFunctional.(c1, c2, 100.0, 50.0) .<= 1 + 1e-9)
    @test any(isapprox.(df.x, 10 * cos(π / 6); atol=1e-8) .&& isapprox.(df.y, 10 * sin(π / 6); atol=1e-8))
    @test !any(isapprox.(df.x, 10.0; atol=1e-8) .&& isapprox.(df.y, 0.0; atol=1e-8))

    # ring
    writeShapeXML(shape_xml, "ellipse", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "100.0", "b" => "50.0", "thickness" => "30.0", "spacing" => string(spacing)]; placement="hexagonal")
    df = generateICCell(shape_xml, shape_domain)
    @test all(ellipseFunctional.(df.x, df.y, 130.0, 80.0) .<= 1 + 1e-9)
    @test all(ellipseFunctional.(df.x, df.y, 100.0, 50.0) .>= 1 - 1e-9)
    @test isapprox(size(df, 1), π * (130 * 80 - 100 * 50) / (spacing^2 * sqrt(3) / 2); rtol=0.03)
end

@testset "rotation parsing" begin
    function parseRotationString(value::Union{Nothing,String}, units::Union{Nothing,String}=nothing)
        e_patch = new_element("patch")
        if !isnothing(value)
            e_rotation = new_child(e_patch, "rotation")
            set_content(e_rotation, value)
            isnothing(units) || set_attribute(e_rotation, "units", units)
        end
        return PhysiCellCellCreator.parseRotation(e_patch, "test patch")
    end
    @test parseRotationString(nothing) == 0.0
    @test parseRotationString("30") ≈ π / 6
    @test parseRotationString("-45", "deg") ≈ -π / 4
    @test parseRotationString("1.5", "rad") ≈ 1.5
    @test parseRotationString("1.5", "radians") ≈ 1.5
    @test parseRotationString("pi") ≈ π
    @test parseRotationString("π") ≈ π
    @test parseRotationString("pi/6") ≈ π / 6
    @test parseRotationString(" π / 6 ") ≈ π / 6
    @test parseRotationString("2pi/3") ≈ 2π / 3
    @test parseRotationString("2*pi/3") ≈ 2π / 3
    @test parseRotationString("-pi/2") ≈ -π / 2
    @test parseRotationString("0.5pi") ≈ π / 2
    @test parseRotationString(".5*π", "rad") ≈ π / 2
    @test_throws AssertionError parseRotationString("pi/6", "deg")
    @test_throws AssertionError parseRotationString("30", "px")
    @test_throws AssertionError parseRotationString("pi^2")
    @test_throws AssertionError parseRotationString("pi/0")
    @test_throws AssertionError parseRotationString("2*3pi")
    @test_throws AssertionError parseRotationString("pi*2")
    @test_throws ArgumentError parseRotationString("thirty")
end

@testset "supported types and elements" begin
    @test PhysiCellCellCreator.supportedPatchTypes() == ["disc", "annulus", "rectangle", "everywhere", "elliptical_disc", "ellipse"]
    @test PhysiCellCellCreator.supportedCarveoutTypes() == ["disc", "annulus", "rectangle", "elliptical_disc", "ellipse"]
    @test PhysiCellCellCreator.supportedPatchTextElements("ellipse") == ["x0", "y0", "z0", "a", "b", "thickness", "normal", "rotation", "number", "spacing", "max_fails", "restrict_to_domain"]
    @test PhysiCellCellCreator.supportedPatchTextElements("elliptical_disc") == ["x0", "y0", "z0", "a", "b", "normal", "rotation", "number", "spacing", "max_fails", "restrict_to_domain"]
    @test "rotation" in PhysiCellCellCreator.supportedPatchTextElements("disc")
    @test !("normal" in PhysiCellCellCreator.supportedPatchTextElements("rectangle"))
    @test PhysiCellCellCreator.supportedPatchTextElements("everywhere") == ["number", "spacing", "max_fails"]
    @test PhysiCellCellCreator.supportedCarveoutTextElements("ellipse") == ["x0", "y0", "z0", "a", "b", "thickness", "normal", "rotation"]
    @test PhysiCellCellCreator.supportedCarveoutTextElements("rectangle") == ["x0", "y0", "z0", "x1", "y1", "width", "height"]
    @test_throws AssertionError PhysiCellCellCreator.supportedPatchTextElements("sphere")
    @test_throws AssertionError PhysiCellCellCreator.supportedCarveoutTextElements("everywhere")
    # unsupported carveout type in an XML
    writeShapeXML(shape_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "10"]; carveouts=[("everywhere", Pair{String,String}[])])
    @test_throws AssertionError generateICCell(shape_xml, shape_domain)
end

@testset "invalid geometry" begin
    for (patch_type, elements) in [
        ("elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "0.0", "b" => "40.0"]),
        ("elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "a" => "40.0", "b" => "-1.0"]),
        ("elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "a" => "40.0", "b" => "40.0"]), # missing z0
        ("elliptical_disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "b" => "40.0"]), # missing a
        ("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "0.0"]),
        ("annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "-1.0", "outer_radius" => "10.0"]),
    ]
        writeShapeXML(shape_xml, patch_type, [elements; "number" => "10"])
        @test_throws AssertionError generateICCell(shape_xml, shape_domain)
    end
end

@testset "template with elliptical patches" begin
    template_folder = joinpath(shape_folder, "template")
    createICCellXMLTemplate(template_folder)
    df = generateICCell(joinpath(template_folder, "cells.xml"), shape_domain)
    # the template's ellipse shell (hexagonal) around the elliptical disc at (-200, 200), both rotated by pi/6
    c1, c2 = rotateBack(df.x, df.y, -200.0, 200.0, π / 6)
    in_shell = ellipseFunctional.(c1, c2, 100.0, 60.0) .<= 1 .&& ellipseFunctional.(c1, c2, 80.0, 40.0) .>= 1
    in_interior = ellipseFunctional.(c1, c2, 80.0, 40.0) .< 1
    @test count(in_shell) > 0
    @test count(in_interior) == 60 # the everywhere patch carves out this ellipse, so only the elliptical_disc patch places cells here
end

rm(shape_folder; recursive=true)
