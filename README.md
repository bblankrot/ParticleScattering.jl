# ParticleScattering

[![Travis](https://travis-ci.org/bblankrot/ParticleScattering.jl.svg?branch=master)](https://travis-ci.org/bblankrot/ParticleScattering.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/p0p636vtrx95ch8m/branch/master?svg=true)](https://ci.appveyor.com/project/bblankrot/particlescattering-jl/branch/master)
[![codecov.io](http://codecov.io/github/bblankrot/ParticleScattering.jl/coverage.svg?branch=master)](http://codecov.io/github/bblankrot/ParticleScattering.jl?branch=master)
[![doc-latest](https://img.shields.io/badge/docs-latest-blue.svg)](https://bblankrot.github.io/ParticleScattering.jl/latest)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00691/status.svg)](https://doi.org/10.21105/joss.00691)

A Julia package for solving large-scale electromagnetic
scattering problems in two dimensions; specifically,
those containing a large number of penetrable smooth
particles. Provides the ability to optimize over the
particle parameters for various design problems.

### Installation

ParticleScattering can be installed using `Pkg.add`. Currently, only Julia 0.6 is supported.

```julia
Pkg.add("ParticleScattering")
using ParticleScattering
```

### Community

The easiest way to contribute is by opening issues! Of course, we'd be more than happy if you implement any fixes and send a PR.
If you have any relevant scattering problems that would make good examples for the docs, feel free to open an issue for that as well.

### Citation

If you publish work that utilizes ParticleScattering, please cite it using:
```
@article{Blankrot2018joss,
  title={ParticleScattering: Solving and optimizing multiple-scattering problems in {Julia}},
  author={Blankrot, Boaz and Heitzinger, Clemens},
  journal={Journal of Open Source Software},
  publisher={The Open Journal},
  volume={3},
  pages={691},
  number={25},
  DOI={10.21105/joss.00691},
  year={2018},
  month={May}
}
```
