# using Pkg
# Pkg.activate(".")


using Documenter, RandomizedProgressiveHedging

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"      # Solve url problem on laptop
    ),
    sitename="RandomizedProgressiveHedging",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "quickstart.md",
        "Library" => Any[
            "Public" => "public_api.md",
            "Internal" => "internal_api.md",
        ],
    ],
    modules = [RandomizedProgressiveHedging]
)

deploydocs(
    repo = "github.com/yassine-laguel/RandomizedProgressiveHedging.jl.git",
)
