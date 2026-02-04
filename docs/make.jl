using Documenter
using MonteCarloGlauber

DocMeta.setdocmeta!(MonteCarloGlauber, :DocTestSetup, :(using MonteCarloGlauber); recursive = true)

makedocs(
    modules = [MonteCarloGlauber],
    sitename = "MonteCarloGlauber.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        collapselevel = 1,
    ),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Quickstart" => "quickstart.md",
        "Examples" => "examples.md",
        "Background" => "background.md",
        "API" => "api.md",
    ],
    doctest = false,
)

repo = get(ENV, "GITHUB_REPOSITORY", "AndreasKirchner/MonteCarloGlauber.jl")
deploydocs(
    repo = "github.com/$(repo).git",
    devbranch = "main",
    push_preview = true,
)
