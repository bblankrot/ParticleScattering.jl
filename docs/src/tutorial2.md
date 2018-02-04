# Tutorial 2: Accelerating Solutions with FMM

In this tutorial, we examine scattering from several hundreds of particles, and
use the built-in Fast Multipole Method (FMM) implementation to provide faster
results.
FMM groups nearby particles and approximates their cumulative effect on "far"
particles. Here, the grouping is done by drawing a square grid over computational
region encompassing the particles. The error of these approximations is fairly
 controllable, especially for mid- to high-frequency scattering.

Another difference between the direct and FMM solvers is that the latter requires
an iterative solver for the system of equations (GMRES is used here). Thus a
certain residual tolerance must be defined at which point the GMRES process
terminates.

The scattering problem here is given by:

```julia
λ0 = 1 #doesn't matter since everything is normalized to λ0
k0 = 2π/λ0
kin = 3k0
θ_i = π/4 #incident wave e^{i k_0 (1/sqrt{2},1/sqrt{2}) \cdot \mathbf{r}}

N_squircle = 200
N_star = 260
P = 10

M = 20
shapes = [rounded_star(0.1λ0, 0.05λ0, 5, N_star);
            squircle(0.15λ0, N_squircle)]
centers =  square_grid(M, 0.4λ0) #MxM grid with distance 0.4λ0
ids = rand(1:2, M^2)
φs = 2π*rand(M^2) #random rotation angles
sp = ScatteringProblem(shapes, ids, centers, φs)
```

To setup FMM, we use the constructor `FMMoptions`, whose options are given by:

```julia
FMM::Bool       # Is FMM used? (default: false)
nx::Integer     # number of groups in x direction (required if dx is not
                # specified)
dx::Real        # group height & width (required if nx is not specified)
acc::Integer    # accuracy digits for translation truncation, and also for
                # GMRES if tol is not given (required)
tol::Real       # GMRES tolerance (default: 10^{-acc})
method::String  # method used: for now only "pre" (default: "pre")
```

For this problem, we choose ``6`` digits of accuracy and grouping into
``(M/2)^2`` boxes. The grouping can be viewed by calling `divideSpace`:

```julia
fmm_options = FMMoptions(true, acc = 6, nx = div(M,2))
divideSpace(centers, fmm_options; drawGroups = true)
```

In this plot, the red markers denote the group centers while stars denote
particle centers (the particles can be drawn on top of this plot with
`draw_shapes`). At first, it might look strange that most each particle lies
outside the FMM group; however, the FMM is used only after the particles are
converted to line sources, and are thus fully contained in the FMM grid.

![fmm_tutorial_plot0](./assets/fmm_tutorial_plot0.png)

Calculating and plotting the near or far fields with FMM is just as in the
[previous tutorial](@ref scattering_small_grid), except we must supply the `FMMoptions` object:

```julia
data = plot_near_field(k0, kin, P, sp, θ_i, opt = fmm_options)
colorbar()
```

![fmm_tutorial_plot1](./assets/fmm_tutorial_plot1.png)
