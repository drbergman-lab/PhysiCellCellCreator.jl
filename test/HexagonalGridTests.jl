using LightXML, DataFrames

hex_folder = "hex_test_folder"
mkpath(hex_folder)

function writePatchXML(path::String, patch_type::String, elements::Vector{Pair{String,String}}; cell_type::String="default", carveout::Union{Nothing,Tuple{String,Vector{Pair{String,String}}}}=nothing)
    xml_doc = XMLDocument()
    xml_root = create_root(xml_doc, "ic_cells")
    e_patches = new_child(xml_root, "cell_patches")
    set_attribute(e_patches, "name", cell_type)
    e_collection = new_child(e_patches, "patch_collection")
    set_attribute(e_collection, "type", patch_type)
    e_patch = new_child(e_collection, "patch")
    set_attribute(e_patch, "ID", "1")
    for (name, value) in elements
        set_content(new_child(e_patch, name), value)
    end
    if !isnothing(carveout)
        carveout_type, carveout_elements = carveout
        e_carveouts = new_child(e_patch, "carveout_patches")
        e_carveout_collection = new_child(e_carveouts, "patch_collection")
        set_attribute(e_carveout_collection, "type", carveout_type)
        e_carveout = new_child(e_carveout_collection, "patch")
        set_attribute(e_carveout, "ID", "1")
        for (name, value) in carveout_elements
            set_content(new_child(e_carveout, name), value)
        end
    end
    save_file(xml_doc, path)
    free(xml_doc)
    return path
end

function nearestNeighborDistances(df::DataFrame)
    n = size(df, 1)
    dists = fill(Inf, n)
    for i in 1:n, j in 1:n
        i == j && continue
        d = sqrt((df.x[i] - df.x[j])^2 + (df.y[i] - df.y[j])^2 + (df.z[i] - df.z[j])^2)
        dists[i] = min(dists[i], d)
    end
    return dists
end

hex_domain = Dict{String,Float64}("x_min" => -500.0, "x_max" => 500.0, "y_min" => -500.0, "y_max" => 500.0, "z_min" => -500.0, "z_max" => 500.0)
hex_xml = joinpath(hex_folder, "cells.xml")

@testset "hexagonal disc" begin
    writePatchXML(hex_xml, "disc", ["x0" => "10.0", "y0" => "-20.0", "z0" => "5.0", "radius" => "40.0", "hexagonal_grid" => "true", "spacing" => "10.0"])
    df = generateICCell(hex_xml, hex_domain)
    @test size(df, 1) > 0
    @test all(df.type .== "default")
    @test all(df.z .== 5.0) # z0 is used exactly once (no doubling)
    r = sqrt.((df.x .- 10.0) .^ 2 .+ (df.y .+ 20.0) .^ 2)
    @test all(r .<= 40.0 + 1e-8)
    @test any(r .< 1e-8) # a cell sits at the center of the disc
    @test all(isapprox.(nearestNeighborDistances(df), 10.0; atol=1e-8))
    # each interior cell has six neighbors at the spacing distance
    center_ind = findfirst(r .< 1e-8)
    d_center = sqrt.((df.x .- df.x[center_ind]) .^ 2 .+ (df.y .- df.y[center_ind]) .^ 2)
    @test count(isapprox.(d_center, 10.0; atol=1e-8)) == 6
end

@testset "hexagonal annulus" begin
    writePatchXML(hex_xml, "annulus", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "inner_radius" => "20.0", "outer_radius" => "60.0", "hexagonal_grid" => "true", "spacing" => "8.0"])
    df = generateICCell(hex_xml, hex_domain)
    @test size(df, 1) > 0
    r = sqrt.(df.x .^ 2 .+ df.y .^ 2)
    @test all(20.0 - 1e-8 .<= r .<= 60.0 + 1e-8)
    @test all(isapprox.(nearestNeighborDistances(df), 8.0; atol=1e-8))
end

@testset "hexagonal disc with normal" begin
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "30.0", "hexagonal_grid" => "true", "spacing" => "10.0", "normal" => "1,0,0"])
    df = generateICCell(hex_xml, hex_domain)
    @test size(df, 1) > 0
    @test all(abs.(df.x) .< 1e-8) # disc lies in the plane x = 0
    @test all(sqrt.(df.y .^ 2 .+ df.z .^ 2) .<= 30.0 + 1e-8)
    @test all(isapprox.(nearestNeighborDistances(df), 10.0; atol=1e-8))
end

@testset "hexagonal rectangle" begin
    writePatchXML(hex_xml, "rectangle", ["x0" => "-50.0", "y0" => "-50.0", "z0" => "0.0", "width" => "100.0", "height" => "100.0", "hexagonal_grid" => "true", "spacing" => "10.0"])
    df = generateICCell(hex_xml, hex_domain)
    @test size(df, 1) > 0
    @test all(-50.0 - 1e-8 .<= df.x .<= 50.0 + 1e-8)
    @test all(-50.0 - 1e-8 .<= df.y .<= 50.0 + 1e-8)
    @test all(isapprox.(nearestNeighborDistances(df), 10.0; atol=1e-8))
    @test any(abs.(df.x) .< 1e-8 .&& abs.(df.y) .< 1e-8) # a cell sits at the center of the rectangle

    # same rectangle with a disc carved out of the center
    writePatchXML(hex_xml, "rectangle", ["x0" => "-50.0", "y0" => "-50.0", "z0" => "0.0", "x1" => "50.0", "y1" => "50.0", "hexagonal_grid" => "true", "spacing" => "10.0"];
        carveout=("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "25.0"]))
    df_carved = generateICCell(hex_xml, hex_domain)
    @test size(df_carved, 1) < size(df, 1)
    @test all(sqrt.(df_carved.x .^ 2 .+ df_carved.y .^ 2) .> 25.0)
end

@testset "hexagonal everywhere" begin
    writePatchXML(hex_xml, "everywhere", ["hexagonal_grid" => "true", "spacing" => "50.0"])
    small_domain = Dict{String,Float64}("x_min" => -100.0, "x_max" => 100.0, "y_min" => -100.0, "y_max" => 100.0, "z_min" => -10.0, "z_max" => 10.0)
    df = generateICCell(hex_xml, small_domain)
    @test size(df, 1) > 0
    @test all(-100.0 - 1e-8 .<= df.x .<= 100.0 + 1e-8)
    @test all(-100.0 - 1e-8 .<= df.y .<= 100.0 + 1e-8)
    @test all(df.z .== 0.0)
    @test all(isapprox.(nearestNeighborDistances(df), 50.0; atol=1e-8))
end

@testset "restrict_to_domain with hexagonal grid" begin
    # disc extends past the domain edge
    small_domain = Dict{String,Float64}("x_min" => 0.0, "x_max" => 100.0, "y_min" => -100.0, "y_max" => 100.0, "z_min" => -10.0, "z_max" => 10.0)
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "30.0", "hexagonal_grid" => "true", "spacing" => "10.0"])
    df_restricted = generateICCell(hex_xml, small_domain)
    @test all(df_restricted.x .>= 0.0)
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "30.0", "hexagonal_grid" => "true", "spacing" => "10.0", "restrict_to_domain" => "false"])
    df_unrestricted = generateICCell(hex_xml, small_domain)
    @test any(df_unrestricted.x .< 0.0)
    @test size(df_unrestricted, 1) > size(df_restricted, 1)
end

@testset "spacing from cell_radius" begin
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "hexagonal_grid" => "true"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain) # no spacing and no cell_radius
    # default spacing_factor is 0.95, so spacing = 0.95 * 2 * cell_radius
    df = generateICCell(hex_xml, hex_domain; cell_radius=6.0)
    @test all(isapprox.(nearestNeighborDistances(df), 0.95 * 12.0; atol=1e-8))
    df = generateICCell(hex_xml, hex_domain; cell_radius=Dict("default" => 4))
    @test all(isapprox.(nearestNeighborDistances(df), 0.95 * 8.0; atol=1e-8))
    @test_throws AssertionError generateICCell(hex_xml, hex_domain; cell_radius=Dict("other_cell_type" => 4.0))
    # explicit spacing_factor
    df = generateICCell(hex_xml, hex_domain; cell_radius=6.0, spacing_factor=1.0)
    @test all(isapprox.(nearestNeighborDistances(df), 12.0; atol=1e-8))
    df = generateICCell(hex_xml, hex_domain; cell_radius=6.0, spacing_factor=0.5)
    @test all(isapprox.(nearestNeighborDistances(df), 6.0; atol=1e-8))
    @test_throws AssertionError generateICCell(hex_xml, hex_domain; cell_radius=6.0, spacing_factor=0.0)

    # explicit spacing takes precedence over cell_radius
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "hexagonal_grid" => "true", "spacing" => "10.0"])
    df = generateICCell(hex_xml, hex_domain; cell_radius=6.0, spacing_factor=1.0)
    @test all(isapprox.(nearestNeighborDistances(df), 10.0; atol=1e-8))

    # cell_radius is ignored by random placement
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "20"])
    df = generateICCell(hex_xml, hex_domain; cell_radius=6.0)
    @test size(df, 1) == 20
end

@testset "incompatible and invalid hexagonal settings" begin
    # number and hexagonal_grid together
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "10", "hexagonal_grid" => "true", "spacing" => "10.0"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain)
    # spacing without hexagonal_grid
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "10", "spacing" => "10.0"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain)
    # hexagonal_grid explicitly false falls back to random placement using number
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "number" => "10", "hexagonal_grid" => "false"])
    @test size(generateICCell(hex_xml, hex_domain), 1) == 10
    # neither number nor hexagonal_grid
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain)
    # non-positive spacing
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "hexagonal_grid" => "true", "spacing" => "0.0"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain)
    # non-boolean hexagonal_grid
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "hexagonal_grid" => "yes", "spacing" => "10.0"])
    @test_throws AssertionError generateICCell(hex_xml, hex_domain)
end

@testset "empty hexagonal grid warns" begin
    # the whole disc is carved out
    writePatchXML(hex_xml, "disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "40.0", "hexagonal_grid" => "true", "spacing" => "10.0"];
        carveout=("disc", ["x0" => "0.0", "y0" => "0.0", "z0" => "0.0", "radius" => "100.0"]))
    df = @test_logs (:warn, r"placed no cells") generateICCell(hex_xml, hex_domain)
    @test size(df, 1) == 0
end

@testset "hexagonal grid in template" begin
    template_folder = joinpath(hex_folder, "template")
    createICCellXMLTemplate(template_folder)
    df = generateICCell(joinpath(template_folder, "cells.xml"), hex_domain)
    # the template's hexagonal disc is centered at (200, 200) with radius 40 and spacing 12
    # (random patches in the template may also place cells in this region, so check the lattice points directly)
    hasCellAt(x, y) = any(isapprox.(df.x, x; atol=1e-8) .&& isapprox.(df.y, y; atol=1e-8))
    @test hasCellAt(200.0, 200.0)
    @test hasCellAt(212.0, 200.0)
    @test hasCellAt(206.0, 200.0 + 12.0 * sqrt(3) / 2)
    @test !hasCellAt(236.0, 200.0 + 2 * 12.0 * sqrt(3) / 2) # outside the disc (radius ~ 42.3 > 40)
    hex_ind = sqrt.((df.x .- 200.0) .^ 2 .+ (df.y .- 200.0) .^ 2) .<= 40.0 + 1e-8
    @test count(hex_ind) >= 37 # 1 + 6 + 12 + 18 lattice points within 3 spacings of the center
end

rm(hex_folder; recursive=true)
