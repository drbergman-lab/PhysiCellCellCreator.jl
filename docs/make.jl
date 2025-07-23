using PhysiCellCellCreator
using Documenter

DocMeta.setdocmeta!(PhysiCellCellCreator, :DocTestSetup, :(using PhysiCellCellCreator); recursive=true)

makedocs(;
    modules=[PhysiCellCellCreator],
    authors="Daniel Bergman <danielrbergman@gmail.com> and contributors",
    sitename="PhysiCellCellCreator.jl",
    format=Documenter.HTML(;
        canonical="https://drbergman-lab.github.io/PhysiCellCellCreator.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Guide" => "man/guide.md"
        ]
    ],
)

deploydocs(;
    repo="github.com/drbergman-lab/PhysiCellCellCreator.jl",
    devbranch="development",
)
