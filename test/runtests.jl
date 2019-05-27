using ParticleScattering
using Test, LinearAlgebra
import Optim, LineSearches

include("scatteredfield_test.jl")
include("multipole_test.jl")
include("fmm_test.jl")
include("optimize_radius_test.jl")
include("optimize_angle_test.jl")
include("minimum_N_P_test.jl")
include("utilities_test.jl")
