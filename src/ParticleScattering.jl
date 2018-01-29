module ParticleScattering

#Core functionality
using IterativeSolvers, LinearMaps, Optim
#For plotting with PyPlot
using PyPlot, PyCall
@pyimport matplotlib.patches as patch #circles,polygons
#For plotting with pgfplots
import DataFrames, CSV, PGFPlotsX; const pgf = PGFPlotsX

include("PS_types.jl")
include("shapes.jl")
include("scattering.jl")
include("multipole.jl")
include("visualization.jl")
include("visualization_pgf.jl")
include("fmm_matrices.jl")
include("fmm_mvp.jl")
include("fmm_main.jl")
include("optimize_phis.jl")
include("optimize_rs.jl")
include("utilities.jl")
include("minimum_N_P.jl")

#methods, shapes.jl
export rounded_star, squircle, square_grid, randpoints, luneburg_grid,
    verify_min_distance
#methods, multipole.jl
export solve_particle_scattering, scattered_field_multipole
#methods, minimum_N_P.jl
export minimumN, minimumP
#methods, visualization.jl
export plot_far_field, plot_near_field, calc_near_field, draw_shapes
#methods, fmm_main.jl
export solve_particle_scattering_FMM

### documented till here####


#methods, scattering.jl
export solvePotentialShape, solvePotentialShapePW,
    scatteredField, solvePotential_forError

#methods, optimize_phis.jl
export optimize_φ_grad
#methods, optimize_rs.jl
export optimize_radius
#methods, visualization_pgf.jl
export plotNearField_pgf, drawShapes_pgf


#types, PS_types.jl
export ScatteringProblem, OptimBuffer, FMMoptions, R_multipole,
    ShapeParams, CircleParams, AbstractShapeParams

#temp
include("optimize_rs_old.jl")
export divideSpace, FMMtruncation, particleExpansion, FMMbuildMatrices

export FMMbuffer, FMMmatrices

end
