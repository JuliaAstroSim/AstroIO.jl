"""
Compile with:
julia --project=docs/ --color=yes docs/make.jl

Generate key:
DocumenterTools.genkeys(user="JuliaAstroSim", repo="git@github.com:JuliaAstroSim/AstroIO.jl.git")
"""

using Documenter

using AstroIO

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

makedocs(
    modules = [AstroIO],
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaastrosim.github.io/AstroIO.jl/dev/",
        assets = ["assets/alpha_small.ico"],
        analytics = "UA-153693590-1",
        highlights = ["llvm", "yaml"],
    ),
    clean = false,
    sitename = "AstroIO.jl",
    authors = "islent",
    # `linkcheck = true` shells out to curl and fails the build when the
    # host cannot reach the internet (e.g. behind a firewall / sandbox).
    # Enable explicitly with `julia docs/make.jl checklinks`, or skip
    # explicitly with `julia docs/make.jl skiplinks`. The default is off
    # so the local build never fails for network reasons.
    linkcheck = "checklinks" in ARGS,
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "manual/guide.md",
            "manual/gadget2.md",
            "manual/csv.md",
            "manual/jld2.md",
            "manual/hdf5.md",
            "manual/houdini.md",
            "manual/confparser.md",
            "manual/tools.md",
        ],
        "Library" => Any[
            "lib/Methods.md",
        ],
    ],
)

deploydocs(
    repo = "github.com/JuliaAstroSim/AstroIO.jl.git",
    target = "build",
)
