using Documenter, CTDE

makedocs(;
    modules=[CTDE],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/adolgert/CTDE.jl/blob/{commit}{path}#L{line}",
    sitename="CTDE.jl",
    authors="Drew Dolgert <adolgert@uw.edu>",
    assets=String[],
)

deploydocs(;
    repo="github.com/adolgert/CTDE.jl",
)
