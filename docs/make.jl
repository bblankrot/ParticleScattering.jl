using Documenter, ParticleScattering

makedocs(format = :html,
    sitename = "ParticleScattering.jl",
    authors = "Boaz Blankrot",
    linkcheck = !("skiplinks" in ARGS),
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial1.md",
            "tutorial2.md"
        ],
        "Choosing Minimal N and P" => "minimalNP.md",
        "api.md"
    ]
)

deploydocs(
    repo   = "github.com/bblankrot/ParticleScattering.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
