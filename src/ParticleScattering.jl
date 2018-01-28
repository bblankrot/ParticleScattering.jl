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
export rounded_star, squircle, square_grid, randpoints,
    luneburg_grid, verify_min_distance
#methods, multipole.jl
export solveParticleScattering, scatteredFieldMultipole
#methods, minimum_N_P.jl
export minimumN, minimumP

### documented till here####

#methods, scattering.jl
export solvePotentialShape, solvePotentialShapePW,
    scatteredField, solvePotential_forError
#methods, visualization.jl
export plotFarField, plotNearField, calculateNearField, drawShapes
#methods, fmm_main.jl
export solveParticleScattering_FMM
#methods, optimize_phis.jl
export optimize_φ_grad
#methods, optimize_rs.jl
export optimize_radius
#methods, visualization_pgf.jl
export plotNearField_pgf, drawShapes_pgf


#types, PS_types.jl
export ScatteringProblem, ShapeParams, CircleParams, AbstractShapeParams,
    FMMoptions, OptimBuffer, R_multipole

#temp
include("optimize_rs_old.jl")
export divideSpace, FMMtruncation, particleExpansion, FMMbuildMatrices

export FMMbuffer, FMMmatrices

end
