using Documenter, InsarTimeseries

makedocs(
    modules = [InsarTimeseries],
    format = Documenter.HTML(),
    checkdocs = :exports,
    sitename = "InsarTimeseries.jl",
    pages = Any["index.md"],
)

deploydocs(repo = "github.com/scottstanie/InsarTimeseries.jl.git")
