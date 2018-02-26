using Documenter, ParticleScattering

makedocs(format = :html,
    sitename = "ParticleScattering.jl",
    authors = "Boaz Blankrot"
    linkcheck = !("skiplinks" in ARGS),
    pages = Any[
        "Main" => "index.md",
        "Tutorials" => Any[
            "tutorial1.md",
            "tutorial2.md"
        ],
        "api.md"
    ]
)

deploydocs(
    repo   = "github.com/bblankrot/ParticleScattering.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)
