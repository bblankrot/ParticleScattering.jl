# ParticleScattering.jl

A Julia package for solving large-scale electromagnetic
scattering problems in two dimensions; specifically,
those containing a large number of penetrable smooth
particles. Provides the ability to optimize over particle parameters for various design problems.

## Installation

ParticleScattering can be installed using `Pkg.add`. Currently, only Julia 0.6 is supported.

```julia
Pkg.add("ParticleScattering")
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
introduction to this package yet cover most of its functionality.
Complex and complete examples involving optimization are available in the
 [examples](https://github.com/bblankrot/ParticleScattering.jl/tree/master/examples)
folder.
