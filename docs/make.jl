using Documenter, ParticleScattering

makedocs(format = :html,
    sitename = "ParticleScattering.jl",
    authors = "Boaz Blankrot",
    linkcheck = !("skiplinks" in ARGS),
    pages = Any[
        "Home" => "index.md",
        "Tutorials" => Any[
            "tutorial1.md",
            "tutorial2.md",
            "tutorial_optim_angle.md",
            "tutorial_optim_radius.md"
        ],
        "Choosing Minimal N and P" => "minimalNP.md",
        "Incident Field Types" => "incident_fields.md",
        "Adding New Shapes" => "new_shapes.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo   = "github.com/bblankrot/ParticleScattering.jl.git",
    target = "build",
    julia = "0.6",
    deps   = nothing,
    make   = nothing
)
