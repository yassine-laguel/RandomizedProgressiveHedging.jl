using Pkg
Pkg.activate(".")


using Documenter, RPH

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"      # Solve url problem on laptop
    ),
    sitename="RPH",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "quickstart.md",
        "Library" => Any[
            "Public" => "public_api.md",
            "Internal" => "internal_api.md",
        ],
    ],
    modules = [RPH]
)

# deploydocs(
#     repo = "github.com/yassine-laguel/RPH.jl.git",
# )