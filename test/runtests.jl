using PhysiCellCellCreator, Test

@testset "PhysiCellCellCreator.jl" begin
    # Write your tests here.
    include("./WriteCellsTests.jl")
    include("./HexagonalGridTests.jl")
    include("./ShapeTests.jl")
end
