using LightXML, DataFrames, LinearAlgebra, Random

Random.seed!(20261015)

line_folder = "line_test_folder"
mkpath(line_folder)
line_xml = joinpath(line_folder, "cells.xml")
line_domain = Dict{String,Float64}("x_min" => -500.0, "x_max" => 500.0, "y_min" => -500.0, "y_max" => 500.0, "z_min" => -500.0, "z_max" => 500.0)

function writeLineXML(path::String, elements::Vector{Pair{String,String}}; patch_type::String="line", placement::Union{Nothing,String}=nothing, carveouts::Vector=[])
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

# parameter t of each cell along the segment p1 -> p2 (0 at p1, 1 at p2) and its distance from the line
function lineParameters(df::DataFrame, p1::Vector{Float64}, p2::Vector{Float64})
    d = p2 .- p1
    L = norm(d)
    t = ((df.x .- p1[1]) .* d[1] .+ (df.y .- p1[2]) .* d[2] .+ (df.z .- p1[3]) .* d[3]) ./ L^2
    off = [norm([df.x[i], df.y[i], df.z[i]] .- (p1 .+ t[i] .* d)) for i in 1:size(df, 1)]
    return t, off
end

function nearestNeighborDistanceLine(df::DataFrame)
    n = size(df, 1)
    dists = fill(Inf, n)
    for i in 1:n, j in 1:n
        i == j && continue
        dists[i] = min(dists[i], sqrt((df.x[i] - df.x[j])^2 + (df.y[i] - df.y[j])^2 + (df.z[i] - df.z[j])^2))
    end
    return dists
end

@testset "line: random placement in 2D" begin
    p1, p2 = [-100.0, 0.0, 0.0], [100.0, 50.0, 0.0]
    writeLineXML(line_xml, ["x1" => "-100.0", "y1" => "0.0", "x2" => "100.0", "y2" => "50.0", "number" => "10000"])
    df = generateICCell(line_xml, line_domain)
    @test size(df, 1) == 10000
    @test all(df.type .== "default")
    @test all(df.z .== 0.0) # z1 and z2 default to 0
    t, off = lineParameters(df, p1, p2)
    @test all(off .< 1e-9)
    @test all(-1e-12 .<= t .<= 1 + 1e-12)
    @test isapprox(count(t .< 0.5) / size(df, 1), 0.5; atol=0.02)
    @test isapprox(count(t .< 0.25) / size(df, 1), 0.25; atol=0.02)
end

@testset "line: random placement in 3D" begin
    p1, p2 = [0.0, 0.0, -100.0], [30.0, 40.0, 100.0]
    writeLineXML(line_xml, ["x1" => "0.0", "y1" => "0.0", "z1" => "-100.0", "x2" => "30.0", "y2" => "40.0", "z2" => "100.0", "number" => "5000"])
    df = generateICCell(line_xml, line_domain)
    t, off = lineParameters(df, p1, p2)
    @test all(off .< 1e-9)
    @test all(-1e-12 .<= t .<= 1 + 1e-12)
    @test minimum(df.z) < -90.0 && maximum(df.z) > 90.0

    # only one z given: the other defaults to 0
    writeLineXML(line_xml, ["x1" => "0.0", "y1" => "0.0", "x2" => "0.0", "y2" => "0.0", "z2" => "100.0", "number" => "1000"])
    df = generateICCell(line_xml, line_domain)
    @test all(abs.(df.x) .< 1e-12) && all(abs.(df.y) .< 1e-12)
    @test all(-1e-12 .<= df.z .<= 100.0 + 1e-12)
    @test minimum(df.z) < 10.0
end

@testset "line: evenly spaced placement" begin
    # along x: cells at -50:10:50, endpoints included
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "spacing" => "10.0"]; placement="hexagonal")
    df = generateICCell(line_xml, line_domain)
    @test size(df, 1) == 11
    @test sort(df.x) ≈ collect(-50.0:10.0:50.0)
    @test all(df.y .== 0.0) && all(df.z .== 0.0)

    # a diagonal line of the same length: same count, same spacing, centered on the midpoint
    p1, p2 = [0.0, 0.0, 0.0], [60.0, 80.0, 0.0]
    writeLineXML(line_xml, ["x1" => "0.0", "y1" => "0.0", "x2" => "60.0", "y2" => "80.0", "spacing" => "10.0"]; placement="hexagonal")
    df = generateICCell(line_xml, line_domain)
    @test size(df, 1) == 11
    t, off = lineParameters(df, p1, p2)
    @test all(off .< 1e-9)
    @test any(isapprox.(t, 0.5; atol=1e-12)) # a cell at the midpoint
    @test any(isapprox.(t, 0.0; atol=1e-12)) && any(isapprox.(t, 1.0; atol=1e-12))
    @test all(isapprox.(nearestNeighborDistanceLine(df), 10.0; atol=1e-8))

    # a spacing that does not divide the length: endpoints are not lattice points
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "spacing" => "15.0"]; placement="hexagonal")
    df = generateICCell(line_xml, line_domain)
    @test sort(df.x) ≈ collect(-45.0:15.0:45.0)

    # a line along z uses the fallback frame
    writeLineXML(line_xml, ["x1" => "5.0", "y1" => "-5.0", "z1" => "-50.0", "x2" => "5.0", "y2" => "-5.0", "z2" => "50.0", "spacing" => "25.0"]; placement="hexagonal")
    df = generateICCell(line_xml, line_domain)
    @test size(df, 1) == 5
    @test all(df.x .== 5.0) && all(df.y .== -5.0)
    @test sort(df.z) ≈ collect(-50.0:25.0:50.0)

    # spacing from cell_radius
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0"]; placement="hexagonal")
    df = generateICCell(line_xml, line_domain; cell_radius=5.0, spacing_factor=1.0)
    @test size(df, 1) == 11
end

@testset "line: carveouts and domain" begin
    # a disc carveout cuts a gap in the line
    carveout = ("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "25.0"])
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "spacing" => "10.0"]; placement="hexagonal", carveouts=[carveout])
    df = generateICCell(line_xml, line_domain)
    @test sort(df.x) ≈ [-50.0, -40.0, -30.0, 30.0, 40.0, 50.0]
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "number" => "1000"]; carveouts=[carveout])
    df = generateICCell(line_xml, line_domain)
    @test size(df, 1) == 1000
    @test all(abs.(df.x) .> 25.0)

    # restrict_to_domain clips the line
    small_domain = Dict{String,Float64}("x_min" => 0.0, "x_max" => 100.0, "y_min" => -100.0, "y_max" => 100.0, "z_min" => -10.0, "z_max" => 10.0)
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "spacing" => "10.0"]; placement="hexagonal")
    df = generateICCell(line_xml, small_domain)
    @test sort(df.x) ≈ collect(0.0:10.0:50.0)
    writeLineXML(line_xml, ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0", "spacing" => "10.0", "restrict_to_domain" => "false"]; placement="hexagonal")
    df = generateICCell(line_xml, small_domain)
    @test size(df, 1) == 11
end

@testset "line: errors and registries" begin
    @test "line" in PhysiCellCellCreator.supportedPatchTypes()
    @test !("line" in PhysiCellCellCreator.supportedCarveoutTypes())
    @test PhysiCellCellCreator.supportedPatchTextElements("line") == ["x1", "y1", "z1", "x2", "y2", "z2", "number", "spacing", "max_fails", "restrict_to_domain"]
    @test_throws AssertionError PhysiCellCellCreator.supportedCarveoutTextElements("line")
    # coincident endpoints
    writeLineXML(line_xml, ["x1" => "10.0", "y1" => "10.0", "x2" => "10.0", "y2" => "10.0", "number" => "10"])
    @test_throws AssertionError generateICCell(line_xml, line_domain)
    # missing endpoint coordinate
    writeLineXML(line_xml, ["x1" => "10.0", "y1" => "10.0", "x2" => "20.0", "number" => "10"])
    @test_throws AssertionError generateICCell(line_xml, line_domain)
    # line as a carveout
    writeLineXML(line_xml, ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "10"]; patch_type="disc",
        carveouts=[("line", ["x1" => "-50.0", "y1" => "0.0", "x2" => "50.0", "y2" => "0.0"])])
    @test_throws AssertionError generateICCell(line_xml, line_domain)
end

@testset "line in template" begin
    template_folder = joinpath(line_folder, "template")
    createICCellXMLTemplate(template_folder)
    df = generateICCell(joinpath(template_folder, "cells.xml"), line_domain)
    # the template's line runs from (-400, -450) to (400, -350) with spacing 20
    t, off = lineParameters(df, [-400.0, -450.0, 0.0], [400.0, -350.0, 0.0])
    on_line = off .< 1e-9 .&& -1e-12 .<= t .<= 1 + 1e-12
    @test count(on_line) == 2 * floor(Int, norm([800.0, 100.0]) / 2 / 20.0) + 1
end

rm(line_folder; recursive=true)
