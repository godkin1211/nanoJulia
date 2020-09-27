using nanoJulia
using Documenter

makedocs(;
    modules=[nanoJulia],
    authors="godkin1211 <jiajyun.ciou@gmail.com> and contributors",
    repo="https://github.com/godkin1211/nanoJulia.jl/blob/{commit}{path}#L{line}",
    sitename="nanoJulia.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://godkin1211.github.io/nanoJulia.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/godkin1211/nanoJulia.jl",
)
