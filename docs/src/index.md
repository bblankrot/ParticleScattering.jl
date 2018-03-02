# ParticleScattering.jl

A Julia package for solving large-scale electromagnetic
scattering problems in two dimensions; specifically,
those containing a large number of penetrable smooth
particles. Provides the ability to optimize over particle parameters for various design problems.

## Installation

Currently, only Julia 0.6 is supported. Once Julia is set up, ParticleScattering can be installed by running
```julia
Pkg.clone(https://github.com/bblankrot/ParticleScattering.jl.git)
using ParticleScattering
```
which also installs the following dependencies:

```julia
IterativeSolvers
LinearMaps
Optim
PyPlot
DataFrames
CSV
PGFPlotsX
```

## Getting Started

Users are encouraged to follow the tutorials, as they provide a gradual
introduction to this package. Complex and complete examples involving
optimization are available in the [examples](https://github.com/bblankrot/ParticleScattering.jl/tree/master/examples)
folder.
