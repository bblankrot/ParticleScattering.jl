var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#ParticleScattering.jl-1",
    "page": "Home",
    "title": "ParticleScattering.jl",
    "category": "section",
    "text": "A Julia package for solving large-scale electromagnetic scattering problems in two dimensions; specifically, those containing a large number of penetrable smooth particles. Provides the ability to optimize over particle parameters for various design problems."
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "ParticleScattering can be installed using Pkg.add. Currently, only Julia 0.6 is supported.Pkg.add(\"ParticleScattering\")\nusing ParticleScatteringwhich also installs the following dependencies:IterativeSolvers\nLinearMaps\nOptim\nPyPlot\nDataFrames\nCSV\nPGFPlotsX"
},

{
    "location": "index.html#Getting-Started-1",
    "page": "Home",
    "title": "Getting Started",
    "category": "section",
    "text": "Users are encouraged to follow the tutorials, as they provide a gradual introduction to this package yet cover most of its functionality. Complex and complete examples involving optimization are available in the  examples folder."
},

{
    "location": "tutorial1.html#",
    "page": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "title": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial1.html#Tutorial-1:-Solving-a-Simple-Multiple-Scattering-Problem-1",
    "page": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "title": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial1.html#Scattering-from-one-particle-1",
    "page": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "title": "Scattering from one particle",
    "category": "section",
    "text": "In order to speed up computations, ParticleScattering depends on grouping identical shapes (rotated versions of the same shape are considered identical). Thus solving multiple-scattering problems is syntactically similar to solving for a single particle.In this example, we simulate TM plane-wave scattering (E_z = e^ik(cos theta_i  sin theta_i) cdot mathbfr) from a single rounded star, parametrized by the equation(x(theta) y(theta)) = R + d cos(5theta)(cos theta sin theta) quad\nR = 01lambda_0 d = 005lambda_0which is supplied by rounded_star. For now, we discretize the shape with N=260 nodes and P=10 cylindrical harmonics &ndash; for more information on the relationship between these parameters and the various resulting errors, see Choosing Minimal N and P.λ0 = 1 #doesn\'t matter since everything is normalized to λ0\nk0 = 2π/λ0\nkin = 3k0\nθ_i = 0.0 #incident wave is left->right\n\nN = 260\nP = 10\nshapes = [rounded_star(0.1λ0, 0.05λ0, 5, N)]\nids = [1] # the particle at centers[1,:] has the parametrization shapes[ids[1]]\ncenters = [0.0 0.0] # our particle is centered at the origin\nφs = [0.0] #zero rotation angle\nsp = ScatteringProblem(shapes, ids, centers, φs)Now that the scattering problem is set up, we solve for the cylindrical harmonics coefficients and potential densities, respectively, by usingbeta,inner = solve_particle_scattering(k0, kin, P, sp::ScatteringProblem, θ_i)These can be used to calculate the scattered field at any point in space using low-level function scatteredfield, or strictly outside the circle of radius shapes[1].R with scattered_field_multipole. For large numbers of calculation points, however, it is easier to use calc_near_field which performs solve_particle_scattering and calculates the total (incident + scattered) field at every point using the most appropriate method:#calculate field on the x-axis passing through the particle\npoints = [linspace(-0.5λ0, 0.5λ0, 200)  zeros(200)]\nu = calc_near_field(k0, kin, P, sp, points, θ_i)We use PyPlot to display the result:using PyPlot\nplot(points[:,1]/λ0, abs.(u))(Image: simple_tutorial_plot1)Similarly, a 2D plot can be drawn of the total field around the scatterer:plot_near_field(k0, kin, P, sp::ScatteringProblem, θ_i;\n    x_points = 201, y_points = 201, border = 0.5λ0*[-1;1;-1;1])(Image: simple_tutorial_plot2)Note: In practice, converting the shape potential densities to cylindrical harmonics is inefficient here as we only have one scatterer, and get_potentialPW would be more accurate. This is meant only as an introductory example to the ParticleScattering syntax."
},

{
    "location": "tutorial1.html#scattering_small_grid-1",
    "page": "Tutorial 1: Solving a Simple Multiple-Scattering Problem",
    "title": "Scattering from a small grid of particles",
    "category": "section",
    "text": "Expanding the example above to a collection of different particles is straightforward:λ0 = 1 #doesn\'t matter since everything is normalized to λ0\nk0 = 2π/λ0\nkin = 3k0\nθ_i = π/2 #incident wave is down->up\n\nN_squircle = 200\nN_star = 260\nP = 10\nshapes = [rounded_star(0.1λ0, 0.05λ0, 5, N_star);\n            squircle(0.15λ0, N_squircle)]\nids = [1;2;2;1]\ncenters =  square_grid(2, 0.4λ0) #2x2 grid with distance 0.4λ0\nφs = 2π*rand(4) #random rotation angles\nsp = ScatteringProblem(shapes, ids, centers, φs)Looking at the 4 times 2 array centers, the coordinates of the m-th shape are given by centers[m,:], and its rotation angle is stored in φs[m]. Likewise, ids[m] tells us if the shape has parametrization shapes[1]  &ndash; in this case a rounded star &ndash; or shapes[2], a squircle. It is  imperative that the order of these arrays remain consistent for the solver to correctly precompute the scattering matrix transformation for each particle. Furthermore, shapes should not contain copies of the same shape, as that will lead to unnecessary computations.Plotting the near field with the codedata = plot_near_field(k0, kin, P, sp, θ_i)\ncolorbar()yields the following near-field plot:(Image: simple_tutorial_plot3)The reason the plot is mostly dark is that plot_near_field automatically scales the colors up to the maximum value calculated, which in this case happens to be an artifact due to inaccurate calculations close to a particle boundary. While this issue can be somewhat alleviated by increasing N, it will remain due to the quadrature method used here. Fortunately, this does not affect the results off the shape boundaries, and can be safely ignored by calling clim([0.0;4.0])."
},

{
    "location": "tutorial2.html#",
    "page": "Tutorial 2: Accelerating Solutions with FMM",
    "title": "Tutorial 2: Accelerating Solutions with FMM",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial2.html#tutorial2-1",
    "page": "Tutorial 2: Accelerating Solutions with FMM",
    "title": "Tutorial 2: Accelerating Solutions with FMM",
    "category": "section",
    "text": "In this tutorial, we examine scattering from several hundreds of particles, and use the built-in Fast Multipole Method (FMM) implementation to provide faster results. FMM groups nearby particles and approximates their cumulative effect on \"far\" particles. Here, the grouping is done by drawing a square grid over computational region encompassing the particles. The error of these approximations is fairly controllable, especially for mid- to high-frequency scattering.Another difference between the direct and FMM solvers is that the latter requires an iterative solver for the system of equations (GMRES is used here). Thus a certain residual tolerance must be defined at which point the GMRES process terminates.The scattering problem here is set up byusing ParticleScattering, PyPlot\nλ0 = 1 #doesn\'t matter since everything is normalized to λ0\nk0 = 2π/λ0\nkin = 3k0\nθ_i = π/4 #incident wave e^{i k_0 (1/sqrt{2},1/sqrt{2}) \\cdot \\mathbf{r}}\n\nN_squircle = 200\nN_star = 210\nP = 10\n\nM = 20\nshapes = [rounded_star(0.1λ0, 0.03λ0, 5, N_star);\n            squircle(0.15λ0, N_squircle)]\ncenters =  square_grid(M, λ0) #MxM grid with distance λ0\nids = rand(1:2, M^2)\nφs = 2π*rand(M^2) #random rotation angles\nsp = ScatteringProblem(shapes, ids, centers, φs)To setup FMM, we use the constructor FMMoptions with the following options:FMM::Bool       # Is FMM used? (default: false)\nnx::Integer     # number of groups in x direction (required if dx is not\n                # specified)\ndx::Real        # group height & width (required if nx is not specified)\nacc::Integer    # accuracy digits for translation truncation, and also for\n                # GMRES if tol is not given (required)\ntol::Real       # GMRES tolerance (default: 10^{-acc})\nmethod::String  # method used: for now only \"pre\" (default: \"pre\")For this problem, we choose 6 digits of accuracy and grouping into (M2)^2 boxes. The grouping can be viewed by calling divideSpace:fmm_options = FMMoptions(true, acc = 6, nx = div(M,2))\ndivideSpace(centers, fmm_options; drawGroups = true)In this plot, the red markers denote the group centers while stars denote particle centers (the particles can be drawn on top of this plot with draw_shapes). At first, it might look strange that parts of many particles lie outside the FMM group; however, the FMM is used only after the particles are converted to line sources, and are thus fully contained in the FMM grid.(Image: fmm_tutorial_plot0)Calculating and plotting the near or far fields with FMM is just as in the previous tutorial, except we must supply the FMMoptions object:plot_near_field(k0, kin, P, sp, θ_i, opt = fmm_options,\n                border = [-12;12;-10;10], x_points = 480, y_points = 400)\ncolorbar()(Image: fmm_tutorial_plot1)(Image: fmm_tutorial_plot2)Note: Currently, FMM is used to accelerate the solution of the scattering problem, but not the field calculation in plot_near_field."
},

{
    "location": "tutorial2.html#Direct-vs.-FMM-timing-1",
    "page": "Tutorial 2: Accelerating Solutions with FMM",
    "title": "Direct vs. FMM timing",
    "category": "section",
    "text": ""
},

{
    "location": "tutorial2.html#Which-is-more-accurate?-1",
    "page": "Tutorial 2: Accelerating Solutions with FMM",
    "title": "Which is more accurate?",
    "category": "section",
    "text": "Intuitively, the direct approach is more accurate than the FMM with its various approximations and iterative solution method. However, inaccuracies can arise in the direct solution of even the simplest scattering problems:k0 = 0.01\nkin = 0.02\nshapes = [squircle(1, 200)]\nids = [1;1]\ncenters = [0.0 0.0; 5.0 0.0]\nphis = [0.0;0.0]\nsp = ScatteringProblem(shapes, ids, centers, phis)\nPmax = 15We solve this problem using the direct approach and with FMM, and then compare both the multipole coefficients beta and the resulting potential densities:betas = Array{Vector}(Pmax)\nbetas_FMM = Array{Vector}(Pmax)\ninners = Array{Vector}(Pmax)\ninners_FMM = Array{Vector}(Pmax)\nfmmopts = ParticleScattering.FMMoptions(true, nx = 1, acc = 9)\nfor P = 1:Pmax\n	betas[P], inners[P] = solve_particle_scattering(k0, kin, P, sp, 0.0;\n                            verbose = false)\n	res, inners_FMM[P] = solve_particle_scattering_FMM(k0, kin, P, sp, 0.0,\n                            fmmopts; verbose = false)\n	betas_FMM[P] = res[1]\nend\n\nerrnorm(x,y) = norm(x-y)/norm(x)\n\nfigure()\nsubplot(2,1,1)\nsemilogy([errnorm(betas_FMM[i], betas[i]) for i = 1:Pmax])\nylabel(\"\\$\\\\Delta \\\\beta\\$\")\nsubplot(2,1,2)\nsemilogy([errnorm(inners_FMM[i], inners[i]) for i = 1:Pmax])\nxlabel(\"\\$ P \\$\")\nylabel(\"\\$ \\\\Delta \\$\" * \" Potential Density\")(Image: direct_vs_fmm0)In both subplots, we see that increasing P actually leads to a decrease in accuracy (plotting the results separately also shows that the FMM results stay the same, while the direct results blow up). This is due to two main reasons - conditioning of the system matrix, and the fact that high-order cylindrical harmonics are responsible for substantially greater potential densities than lower-order ones. Both of these are impacted by the number of particles as well as the wavelength.This ties in with Choosing Minimal N and P &ndash; not only does increasing P far beyond that required for a certain error impact runtime, but can also increase the error in the solution.Of course, FMM was not really used here as nx == 1 means both particles are in the same FMM group, and the maintained accuracy is purely due to the iterative system matrix solution used in solve_particle_scattering_FMM, GMRES."
},

{
    "location": "tutorial_optim_angle.html#",
    "page": "Tutorial 3: Angle Optimization",
    "title": "Tutorial 3: Angle Optimization",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial_optim_angle.html#Tutorial-3:-Angle-Optimization-1",
    "page": "Tutorial 3: Angle Optimization",
    "title": "Tutorial 3: Angle Optimization",
    "category": "section",
    "text": "In this tutorial, we build upon the [previous tutorial]@(ref tutorial2) by optimizing the rotation angles of the particles (φs) to maximize the field intensity at a specific point. Depending on the scattering problem, wavelengths, and incident field, optimization can have a major or minor impact on the field. The simplest way to perform this optimization is by calling optimize_φ, which in turn utilizes Optim, a Julia package for nonlinear optimization. The type of objective function handled here is given byf_mathrmobj = sum_mathbfr in I  u(mathbfr)^2where I is a set of points that lie outside all scattering discs and u is the z-component of the electric field. This function can be minimized directly, or maximized by minimizing -f_mathrmobj.First we set up our scattering problem:λ0 = 1 #doesn\'t matter since everything is normalized to λ0\nk0 = 2π/λ0\nkin = 0.5k0\nθ_i = 0 #incident wave e^{i k_0 (1/sqrt{2},1/sqrt{2}) \\cdot \\mathbf{r}}\nM = 20\nshapes = [rounded_star(0.35λ0, 0.1λ0, 4, 202)]\nP = 12\ncenters =  rect_grid(2, div(M,2), λ0, λ0) #2xM/2 grid\nids = ones(Int, M)\nφs0 = zeros(M)\nsp = ScatteringProblem(shapes, ids, centers, φs0)\nfmm_options = FMMoptions(true, acc = 6, dx = 2λ0)\npoints = 0.05*λ0*[-1 0; 1 0; 0 1; 0 -1]where φs0 is the starting point for the optimization method, and points are the locations at which we intend to maximize or minimize the field intensity. In this case, we want to optimize intensity at a small area around the origin. We now select the optimization method and select its options. In most cases, this combination of BFGS with a backtracking line search will yield accurate results in fast time; other line searches that require re-evaluation of the gradient will be significantly slower but may converge more accurately. The possible convergence criteria are set by the bounds f_tol, g_tol, and x_tol, for a relative change in the function, gradient norm, or variables, respectively. In addition, we can set a maximum number of iterations. Verbosity of the output is set with show_trace and extended_trace.optim_options = Optim.Options(f_tol = 1e-6, iterations = 100,\n                    store_trace = true, show_trace = true)\noptim_method = Optim.BFGS(;linesearch = LineSearches.BackTracking())We now run both minimization and maximization:res_min = optimize_φ(φs0, points, P, θ_i, k0, kin, shapes, centers, ids,\n            fmm_options, optim_options, optim_method; minimize = true)\nres_max = optimize_φ(φs0, points, P, θ_i, k0, kin, shapes, centers, ids,\n            fmm_options, optim_options, optim_method; minimize = false)\nsp_min = ScatteringProblem(shapes, ids, centers, res_min.minimizer)\nsp_max = ScatteringProblem(shapes, ids, centers, res_max.minimizer)Once the optimization is done, we can visualize each ScatteringProblem separately with plot_near_field or compare them side by side with the following PyPlot code:plts = Array{Any}(3)\nplts[1] = plot_near_field(k0, kin, P, sp, θ_i, x_points = 100, y_points = 300,\n        opt = fmm_options, border = find_border(sp, points))\nplts[2] = plot_near_field(k0, kin, P, sp_min, θ_i, x_points = 100, y_points = 300,\n        opt = fmm_options, border = find_border(sp, points))\nplts[3] = plot_near_field(k0, kin, P, sp_max, θ_i, x_points = 100, y_points = 300,\n        opt = fmm_options, border = find_border(sp, points))\nclose(\"all\")\n\nfig, axs = subplots(ncols=3); msh = 0\nfor (i, spi) in enumerate([sp;sp_min;sp_max])\n    msh = axs[i][:pcolormesh](plts[i][2][1], plts[i][2][2], abs.(plts[i][2][3]),\n                        vmin = 0, vmax = 3.4, cmap=\"viridis\")\n    draw_shapes(spi.shapes, spi.centers, spi.ids, spi.φs, axs[i])\n    axs[i][:set_aspect](\"equal\", adjustable = \"box\")\n    axs[i][:set_xlim]([border[1];border[2]])\n    axs[i][:set_ylim]([border[3];border[4]])\n    axs[i][:set_xticks]([-1,0,1])\n    i > 1 && axs[i][:set_yticks]([])\nend\nsubplots_adjust(left=0.05, right=0.8, top=0.98, bottom = 0.05, wspace = 0.1)\ncbar_ax = fig[:add_axes]([0.85, 0.05, 0.05, 0.93])\nfig[:colorbar](msh, cax=cbar_ax)<div style=\"text-align:center\">\n<img alt=optim_angle src=\"./assets/optim_angle.png\" style=\"width:80%; height:auto; margin:1%; max-width: 600px\">\n</div><p style=\"clear:both;\">From left to right, we see the electric field before optimization, after minimization, and after maximization. The field intensity at the origin is notably different in both optimization results, with minimization decreasing the intensity by 95%, and maximization increasing it by over 700%. The convergence of the optimization method for both examples can be plotted via:iters = length(res_min.trace)\nfobj = [res_min.trace[i].value for i=1:iters]\ngobj = [res_min.trace[i].g_norm for i=1:iters]\niters2 = length(res_max.trace)\nfobj2 = -[res_max.trace[i].value for i=1:iters2]\ngobj2 = [res_max.trace[i].g_norm for i=1:iters2]\n\nfig, axs = subplots(ncols=2, figsize=[7,5])\naxs[1][:semilogy](0:iters-1, fobj, linewidth=2)\naxs[2][:semilogy](0:iters-1, gobj, linewidth=2)\naxs[1][:semilogy](0:iters2-1, fobj2, linewidth=2, \"--\")\naxs[2][:semilogy](0:iters2-1, gobj2, linewidth=2, \"--\")\naxs[1][:legend]([\"\\$f_\\\\mathrm{obj}\\$ (min)\";\n                \"\\$f_\\\\mathrm{obj}\\$ (max)\"], loc=\"right\")\naxs[2][:legend]([\"\\$\\\\|\\\\mathbf{g}_\\\\mathrm{obj}\\\\|\\$ (min)\";\n                \"\\$\\\\|\\\\mathbf{g}_\\\\mathrm{obj}\\\\|\\$ (max)\"], loc=\"best\")\naxs[1][:set_xlabel](\"Iteration\")\naxs[2][:set_xlabel](\"Iteration\")\naxs[1][:set_ylim](ymax=40)<p style=\"text-align:center;\"><img alt=optim_angle_conv src=\"./assets/optim_angle_conv.png\" style=\"width:70%; height:auto; max-width:600px\"></p>"
},

{
    "location": "tutorial_optim_radius.html#",
    "page": "Tutorial 4: Radius Optimization",
    "title": "Tutorial 4: Radius Optimization",
    "category": "page",
    "text": ""
},

{
    "location": "tutorial_optim_radius.html#Tutorial-4:-Radius-Optimization-1",
    "page": "Tutorial 4: Radius Optimization",
    "title": "Tutorial 4: Radius Optimization",
    "category": "section",
    "text": "In this tutorial, we explore radius optimization, whereby the radii of a group of circular particles are optimized simultaneously to minimize or maximize the electric field intensity at a group of points.We begin by defining a standard scattering scenario with a 5 times 5 square grid of particles:using PyPlot, ParticleScattering\nimport Optim\n\ner = 4.5\nk0 = 2π\nkin = sqrt(er)*k0\na = 0.2*2π/k0     #wavelength/5\nθ_i = 0.0\nP = 5\ncenters = square_grid(5, a)\nφs = zeros(size(centers,1))\nfmm_options = FMMoptions(true, acc = 6, dx = 2a)optimize_radius not only allows us to optimize all of the radii simultaneously, but also to assign several particles the same id, which can be useful when the target radii are expected to have symmetry of some type. Here we shall assume symmetry with respect to the x-axis (horizontal line of symmetry) with uniqueind:# let\'s impose symmetry wrt x-axis\ncenters_abs = centers[:,1] + 1im*abs.(centers[:,2])\nids, centers_abs = uniqueind(centers_abs)\nJ = maximum(ids) #number of optim varsThe same could be done for the y-axis, both axes simultaneously, or radial symmetry, by appropriately choosing center_abs. We now define the optimization parameters via Optim.Options, with convergence decided by the radii and a limited number of 5 outer iterations (with up to 5 inner iterations each). We choose to minimize the field intensity at a single point outside the structure, assert that this point will remain outside the particles regardless of their size, and set the lower and upper bounds for each circle:optim_options =  Optim.Options(x_tol = 1e-6, iterations = 5,\n                               store_trace = true, show_trace = true,\n                               allow_f_increases = true)\n\npoints = [4a 0.0]\nr_max = (0.4*a)*ones(J)\nr_min = (1e-3*a)*ones(J)\nrs0 = (0.25*a)*ones(J)\nassert(verify_min_distance([CircleParams(r_max[i]) for i = 1:J],\n        centers, ids, points))The optimization process is initiated by running:res = optimize_radius(rs0, r_min, r_max, points, ids, P, θ_i, k0, kin,\n                centers, fmm_options, optim_options, minimize = true)\nrs = res.minimizerWith the optimization process complete, we can plot the electric field with the initial and optimized radii:sp1 = ScatteringProblem([CircleParams(rs0[i]) for i = 1:J], ids, centers, φs)\nplot_near_field(k0, kin, P, sp1, θ_i, x_points = 150, y_points = 150,\n        opt = fmm_options, border = 0.9*[-1;1;-1;1], normalize = a)\ncolorbar()\nclim([0;2.5])\nxlabel(\"x/a\")\nylabel(\"y/a\")\nsp2 = ScatteringProblem([CircleParams(rs[i]) for i = 1:J], ids, centers, φs)\nplot_near_field(k0, kin, P, sp2, θ_i, x_points = 150, y_points = 150,\n        opt = fmm_options, border = 0.9*[-1;1;-1;1], normalize = a)\ncolorbar()\nclim([0;2.5])\nxlabel(\"x/a\")\nylabel(\"y/a\")<div style=\"text-align:center\">\n<img alt=optim_radius_before src=\"./assets/optim_radius_before.png\" style=\"width:40%; height:auto; margin:1%; max-width: 300px\">\n<img alt=optim_radius_after src=\"./assets/optim_radius_after.png\" style=\"width:40%; height:auto; margin:1%; max-width: 300px\">\n</div><p style=\"clear:both;\">res also stores the objective value as well as the g radient norm in each iteration. This can be extracted byinner_iters = length(res.trace)\niters = [res.trace[i].iteration for i=1:inner_iters]\nfobj = [res.trace[i].value for i=1:inner_iters]\ngobj = [res.trace[i].g_norm for i=1:inner_iters]\nrng = iters .== 0where rng now contains the indices at which a new outer iteration has begun. Finally, plotting fobj and gobj for this example yields the following plot:<p style=\"text-align:center;\"><img alt=optim_radius_conv src=\"./assets/optim_radius_conv.png\" style=\"width:60%; height:auto; max-width:400px\"></p>where markers denote the start of an outer iteration."
},

{
    "location": "minimalNP.html#",
    "page": "Choosing Minimal N and P",
    "title": "Choosing Minimal N and P",
    "category": "page",
    "text": ""
},

{
    "location": "minimalNP.html#minimalNP-1",
    "page": "Choosing Minimal N and P",
    "title": "Choosing Minimal N and P",
    "category": "section",
    "text": ""
},

{
    "location": "minimalNP.html#ParticleScattering.minimumN",
    "page": "Choosing Minimal N and P",
    "title": "ParticleScattering.minimumN",
    "category": "function",
    "text": "minimumN(kout, kin, shape_function; tol = 1e-9, N_points = 10_000,\n    N_start = 400, N_min = 100, N_max = 1_000) -> N, err\n\nReturn the minimum N necessary (i.e. 2N nodes) to achieve error of at most tol in the electric field for a ShapeParams inclusion created by shape_function(N) which is filled with material of wavenumber kin and surrounded by free space with wavenumber k0. Error is calculated on N_points points on the scattering disk (s.R), by assuming a fictitious line source and comparing its field to that produced by the resulting potential densities.\n\nSince the error scales with N^-3 for moderate wavelengths and errors, we estimate N using the error of N_start, then binary search based on that guess between N_min and N_max.\n\n\n\n"
},

{
    "location": "minimalNP.html#Discretization-Parameter-N-1",
    "page": "Choosing Minimal N and P",
    "title": "Discretization Parameter N",
    "category": "section",
    "text": "For each non-circular shape in a given scattering problem, we must choose the number of discretization nodes 2N that not only fulfills some accuracy requirement, but also is not large enough to slows down the solution process. Although each shape is only solved or once in the pre-processing stage, with O(N^3) time complexity this stage can be slower than the system matrix solution for large values of N.As the relationship between N and the resulting error depends not only on the geometry and diameter of the shape, but also on the wavelengths inside and outside of it, a general approach to computing N is crucial for dependable results. Moreover, there are many ways to quantify the error for a given discretization. Here we utilize a fictitious-source approach: the potential densities sigma, mu are solved for assuming (for example) a plane wave outside the particle and a line source inside it. The fields induced by these densities outside the particle should then be equal to that of the line source, up to some error. Specifically, this error is measured on a circle of radiusR_multipole*maximum(hypot.(s.ft[:,1],s.ft[:,2]))as this is the scattering disc on which the translation to cylindrical harmonics (Bessel & Hankel functions) are performed (and beyond which any gain or loss of accuracy due to N is mostly irrelevant).Above a certain N, this error tends to decay as O(N^-3), but with a multiplicative factor that is heavily dependent on the particle and wavelength. With minimumN, we first guess a value and then use a binary search to find the minimal N satisfying some error tolerance tol:minimumN"
},

{
    "location": "minimalNP.html#Multipole-Parameter-P-1",
    "page": "Choosing Minimal N and P",
    "title": "Multipole Parameter P",
    "category": "section",
    "text": ""
},

{
    "location": "new_shapes.html#",
    "page": "Adding New Shapes",
    "title": "Adding New Shapes",
    "category": "page",
    "text": ""
},

{
    "location": "new_shapes.html#Adding-New-Shapes-1",
    "page": "Adding New Shapes",
    "title": "Adding New Shapes",
    "category": "section",
    "text": "ParticleScattering includes functions for drawing squircles, rounded stars, and ellipses. New shape functions can  be added, provided they have the following structure:function my_shape(args, N)\n    t = Float64[π*j/N for j = 0:(2*N-1)] # or t = 0:π/N:π*(2-1/N)\n    ft =  [x    y]\n    dft = [dx/dt    dy/dt]\n    ShapeParams(t, ft, dft)\nendWhere t is the parametrization variable, ft[i,:] = [x(t[i]) y(t[i])] contains the coordinates, and dft contains the derivative of ft with respect to t. In particular, the quadrature used by ParticleScattering assumes t are equidistantly distributed in 0 2pi), and that none of the points ft lie on the origin."
},

{
    "location": "api.html#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api.html#Index-1",
    "page": "API",
    "title": "Index",
    "category": "section",
    "text": ""
},

{
    "location": "api.html#ParticleScattering",
    "page": "API",
    "title": "ParticleScattering",
    "category": "module",
    "text": "A Julia package for solving large-scale electromagnetic scattering problems in two dimensions; specifically, those containing a large number of penetrable smooth particles. Provides the ability to optimize over the particle parameters for various design problems.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.R_multipole",
    "page": "API",
    "title": "ParticleScattering.R_multipole",
    "category": "constant",
    "text": "R_multipole = 1.1 is the constant ratio between scattering disks and the maximal radius of their particles, and thus half the minimal distance between neighboring particles. While mathematically this can be reduced to 1 + eps(), that will increase the necessary P.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.AbstractShapeParams",
    "page": "API",
    "title": "ParticleScattering.AbstractShapeParams",
    "category": "type",
    "text": "AbstractShapeParams\n\nAbstract type which all shape types inherit from.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.CircleParams",
    "page": "API",
    "title": "ParticleScattering.CircleParams",
    "category": "type",
    "text": "CircleParams(R)\n\nReturns object for a circular shape, containing its radius in the field R (which is also the radius of the scattering disk).\n\nSee also: ShapeParams,R_multipole.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.FMMoptions",
    "page": "API",
    "title": "ParticleScattering.FMMoptions",
    "category": "type",
    "text": "FMMoptions(FMM; nx = 0, dx = 0.0, acc = 0, tol = 0.0, method = \"pre\")\n\nConstructor for struct containing all FMM options. FMM decides if FMM is used, and the following keyword arguments dictate its behavior:\n\nnx::Integer: number of groups in x direction (for division)\ndx::Real: group height/width (alternative division)\nacc::Integer: accuracy digits for translation truncation, and also for gmres if tol is not given\ntol::Real: gmres tolerance\nmethod::String: method used: for now can be \"pre\" or \"pre2\". Mainly used for development.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.OptimBuffer",
    "page": "API",
    "title": "ParticleScattering.OptimBuffer",
    "category": "type",
    "text": "OptimBuffer(Ns::Integer, P::Integer, Npoints::Integer, [J::Integer])\n\nConstructor for the OptimBuffer type, which stores some of the buffers and shared variables necessary for optimization. Includes the cylindrical harmonics coefficient vector β, field values at points of interest (f), the partial derivatives ∂β, and storage for the various right-hand side vectors used while solving for ∂β.\n\nIf the number of optimization variables J is not supplied, it is assumed that J = Ns.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.ScatteringProblem",
    "page": "API",
    "title": "ParticleScattering.ScatteringProblem",
    "category": "type",
    "text": "ScatteringProblem(shapes, ids, centers, φs)\n\nConstructor for the ScatteringProblem type, including particle shape information for multiple-scattering problems.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.ShapeParams",
    "page": "API",
    "title": "ParticleScattering.ShapeParams",
    "category": "type",
    "text": "ShapeParams(t,ft,dft)\n\nReturns ShapeParams object containing the parametrization of a two-dimensional shape. t is a uniform sampling of [0,2π), ft = [x(t) y(t)], and dft = [x\'(t) y\'(t)]. The field R contains the radius of the shape\'s scattering disk.\n\nSee also: CircleParams,R_multipole.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.calc_near_field-Tuple{Any,Any,Any,ParticleScattering.ScatteringProblem,Any,Any}",
    "page": "API",
    "title": "ParticleScattering.calc_near_field",
    "category": "method",
    "text": "calc_near_field(k0, kin, P, sp::ScatteringProblem, points, θ_i;\n                        opt::FMMoptions = FMMoptions(), use_multipole = true,\n                        verbose = true)\n\nCalculates the total electric field as a result of a plane wave with incident angle θ_i scattering from the ScatteringProblem sp, at points. Uses the FMM options given by opt (default behavious is disabled FMM); use_multipole dictates whether electric field is calculated using the multipole/cylindrical harmonics (true) or falls back on potential densities (false). Either way, the multiple-scattering system is solved in the cylindrical harmonics space, and the field by a particular scatterer inside its own scattering discs is calculated by potential densities, as the cylindrical harmonics approximation is not valid there.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.draw_shapes",
    "page": "API",
    "title": "ParticleScattering.draw_shapes",
    "category": "function",
    "text": "draw_shapes(shapes, centers, ids, φs, ax = gca())\n\nDraws all of the shapes in a given scattering problem. Parametrized shapes are drawn as polygons while circles are drawn using matplotlib\'s patch.Circle.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.ellipse-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "ParticleScattering.ellipse",
    "category": "method",
    "text": "ellipse(r1, r2, N)\n\nReturn a ShapeParams object containing the shape parametrized by (x/r1)^2 + (y/r2)^2 = 1 with 2N nodes.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.find_border-Tuple{ParticleScattering.ScatteringProblem,Array{Float64,2}}",
    "page": "API",
    "title": "ParticleScattering.find_border",
    "category": "method",
    "text": "find_border(sp::ScatteringProblem, points::Array{Float64,2}) -> [x_min; x_max; y_min; y_max]\n\nReturns bounding box that contains all of the shapes in sp as well as specified points.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.find_border-Tuple{ParticleScattering.ScatteringProblem}",
    "page": "API",
    "title": "ParticleScattering.find_border",
    "category": "method",
    "text": "find_border(sp::ScatteringProblem) -> [x_min; x_max; y_min; y_max]\n\nReturns bounding box that contains all of the shapes in sp.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.get_potential-NTuple{6,Any}",
    "page": "API",
    "title": "ParticleScattering.get_potential",
    "category": "method",
    "text": "get_potential(kout, kin, P, s::ShapeParams) -> sigma_mu\n\nGiven a shape s with 2N discretization nodes, outer and inner wavenumbers kout,kin, and the cylindrical harmonics parameter P, returns the potential densities sigma_mu. Each column contains the response to a different harmonic, where the first 2N entries contain the single-layer potential density (sigma), and the lower entries contain the double-layer density (mu).\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.get_potential-Tuple{Any,Any,Any,ParticleScattering.ShapeParams}",
    "page": "API",
    "title": "ParticleScattering.get_potential",
    "category": "method",
    "text": "get_potential(kout, kin, P, t, ft, dft) -> sigma_mu\n\nSame, but with the ShapeParams supplied directly.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.get_potentialPW-NTuple{4,Any}",
    "page": "API",
    "title": "ParticleScattering.get_potentialPW",
    "category": "method",
    "text": "get_potentialPW(kout, kin, s::ShapeParams, θ_i) -> sigma_mu\n\nGiven a shape s with 2N discretization nodes, outer and inner wavenumbers kout,kin, and an incident plane-wave angle, returns the potential densities vector sigma_mu. The first 2N entries contain the single-layer potential density (sigma), and the lower entries contain the double-layer density (mu).\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.hex_grid-Tuple{Integer,Integer,Any}",
    "page": "API",
    "title": "ParticleScattering.hex_grid",
    "category": "method",
    "text": "hex_grid(a::Integer, b::Integer, d)\n\nReturn centers, an (M,2) array  containing the points on a hexagonal lattice with horizontal rows, with a points in each row and rows rows, distanced d. If minus1 is true, the last point in every odd row is omitted.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.luneburg_grid-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "ParticleScattering.luneburg_grid",
    "category": "method",
    "text": "luneburg_grid(R_lens, N_cells, er; levels = 0, TM = true) -> centers, ids, rs\n\nReturns the coordinates and radii of the circular inclusions in a Luneburg lens device of radius R_lens with N_cells unit cells across its diameter. Radii are determined by averaging over cell permittivity, assuming air outside and relative permittivity er in the rods, and depends on incident field polarization (TM/TE with respect to z-axis). If levels == 0, groups identical radii together, such that rs[ids[n]] is the radius of the rod centered at (center[n,1],center[n,2]). Otherwise quantizes the radii to uniformly spaced levels.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.minimumN-Tuple{Any,Any,Any}",
    "page": "API",
    "title": "ParticleScattering.minimumN",
    "category": "method",
    "text": "minimumN(kout, kin, shape_function; tol = 1e-9, N_points = 10_000,\n    N_start = 400, N_min = 100, N_max = 1_000) -> N, err\n\nReturn the minimum N necessary (i.e. 2N nodes) to achieve error of at most tol in the electric field for a ShapeParams inclusion created by shape_function(N) which is filled with material of wavenumber kin and surrounded by free space with wavenumber k0. Error is calculated on N_points points on the scattering disk (s.R), by assuming a fictitious line source and comparing its field to that produced by the resulting potential densities.\n\nSince the error scales with N^-3 for moderate wavelengths and errors, we estimate N using the error of N_start, then binary search based on that guess between N_min and N_max.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.minimumP-Tuple{Any,Any,ParticleScattering.ShapeParams}",
    "page": "API",
    "title": "ParticleScattering.minimumP",
    "category": "method",
    "text": "minimumP(k0, kin, s::ShapeParams; tol = 1e-9, N_points = 10_000, P_min = 1,\n    P_max = 60, dist = 2) -> P, errP\n\nReturn the minimum P necessary to achieve error of at most tol in the electric field, when compared to that obtained with 2N discretization, for a ShapeParams inclusion filled with material of wavenumber kin and surrounded by free space with wavenumber k0. Error is calculated on N_points points on a disk of radius dist*s.R.\n\nUses binary search between P_min and P_max.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.optimize_radius-Tuple{Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Optim.Options}",
    "page": "API",
    "title": "ParticleScattering.optimize_radius",
    "category": "method",
    "text": "optimize_radius(rs0, r_min, r_max, points, ids, P, θ_i, k0, kin, centers,\n    fmmopts, optimopts::Optim.Options; minimize = true, method = \"BFGS\")\n\nOptimize the radii of circular particles for minimization or maximization of the field intensity at points, depending on minimize. Uses Optim\'s Fminbox box-contrained optimization to contain radii in feasible rangle, given in scalar or vector form by r_min and r_max.\n\nHere, ids allows for grouping particles - for example, to maintain symmetry of the optimized device. optimopts defines the convergence criteria and other optimization parameters for both the inner and outer iterations. method can be either \"BFGS\" or \"LBFGS\". See the Optim.Fminbox documentation for more details.\n\nReturns an object of type Optim.MultivariateOptimizationResults.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.optimize_φ-Tuple{Any,Any,Any,Any,Any,Any,Any,Any,Any,Any,Optim.Options,Any}",
    "page": "API",
    "title": "ParticleScattering.optimize_φ",
    "category": "method",
    "text": "optimize_φ(φs0, points, P, θ_i, k0, kin, shapes, centers, ids, fmmopts,\n    optimopts::Optim.Options, minimize = true)\n\nOptimize the rotation angles of a particle collection for minimization or maximization (depending on minimize) of the field intensity at points. optimopts and method define the optimization method, convergence criteria, and other optimization parameters. Returns an object of type Optim.MultivariateOptimizationResults.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.plot_far_field",
    "page": "API",
    "title": "ParticleScattering.plot_far_field",
    "category": "function",
    "text": "plot_far_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;\n                    opt::FMMoptions = FMMoptions(), use_multipole = true,\n                    plot_points = 200)\n\nPlots the echo width (radar cross section in two dimensions) for a given scattering problem. opt, use_multipole are as in plot_near_field. Also returns the echo width.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.plot_near_field",
    "page": "API",
    "title": "ParticleScattering.plot_near_field",
    "category": "function",
    "text": "plot_near_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;\n                    opt::FMMoptions = FMMoptions(), use_multipole = true,\n                    x_points = 201, y_points = 201, border = find_border(sp),\n                    normalize = 1.0)\n\nPlots the total electric field as a result of a plane wave with incident angle θ_i scattering from the ScatteringProblem sp, using matplotlib\'s pcolormesh. Can accept number of sampling points in each direction plus bounding box or calculate automatically.\n\nUses the FMM options given by opt (FMM is disabled by default); use_multipole dictates whether electric field is calculated using the multipole/cylindrical harmonics (true) or falls back on potential densities (false). Either way, the multiple-scattering system is solved in the cylindrical harmonics space. Normalizes all distances and sizes in plot (but not output) by normalize.\n\nReturns the calculated field in two formats:\n\n(points, Ez) where Ez[i] is the total electric field at points[i,:], and\n(xgrid,ygrid,zgrid), the format suitable for pcolormesh, where zgrid[i,j]\n\ncontains the field at (mean(xgrid[i, j:j+1]), mean(ygrid[i:i+1, j])).\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.plot_near_field_pgf",
    "page": "API",
    "title": "ParticleScattering.plot_near_field_pgf",
    "category": "function",
    "text": "plot_near_field_pgf(filename, k0, kin, P, sp::ScatteringProblem, θ_i = 0;\n                    opt::FMMoptions = FMMoptions(), use_multipole = true,\n                    x_points = 201, y_points = 201, border = find_border(sp),\n                    downsample = 1, include_preamble = false, normalize = 1.0)\n\nPlots the total electric field as a result of a plane wave with incident angle θ_i scattering from the ScatteringProblem sp, using pgfplots\'s surf. Can accept number of sampling points in each direction, and either a given border or calculate it automatically. The plots of the shapes (but not the field) can be downsampled by setting an integer downsample, since pgfplots slows down dramatically when drawing many shapes with many nodes.\n\nUses the FMM options given by opt (FMM is disabled by default); use_multipole dictates whether electric field is calculated using the multipole/cylindrical harmonics (true) or falls back on potential densities (false). Either way, the multiple-scattering system is solved in the cylindrical harmonics space. Normalizes all distances and sizes in plot by normalize.\n\nSaves the generated pgfplots file to filename, with just a surrounding tikzpicture environment if include_preamble=false, and a compilable tandalone document otherwise.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.randpoints-NTuple{4,Any}",
    "page": "API",
    "title": "ParticleScattering.randpoints",
    "category": "method",
    "text": "randpoints(M, dmin, width, height; failures = 100)\n\nReturn centers, an (M,2) array containing M points distanced at least dmin in a width by height box. Fails failures times successively before giving up.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.randpoints-NTuple{5,Any}",
    "page": "API",
    "title": "ParticleScattering.randpoints",
    "category": "method",
    "text": "randpoints(M, dmin, width, height, points; failures = 100)\n\nSame as randpoints(M, dmin, width, height; failures = 100) but also requires centers to be distanced at least dmin from points.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.rect_grid-Tuple{Integer,Integer,Any,Any}",
    "page": "API",
    "title": "ParticleScattering.rect_grid",
    "category": "method",
    "text": "rect_grid(a::Integer, b::Integer, dx, dy)\n\nReturn centers, an (a*b,2) array containing the points spanned by a points distanced dx and b points distanced dy, in the x and y directions, respectively.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.rounded_star-NTuple{4,Any}",
    "page": "API",
    "title": "ParticleScattering.rounded_star",
    "category": "method",
    "text": "rounded_star(r, d, num, N)\n\nReturn a ShapeParams object containing the shape parametrized by (x()y()) = (r + d*cos(*num))*(cos()sin()) with 2N nodes.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.scatteredfield-NTuple{6,Any}",
    "page": "API",
    "title": "ParticleScattering.scatteredfield",
    "category": "method",
    "text": "scatteredfield(sigma_mu, k, t, ft, dft, p) -> u_s\n\nSame, but with the ShapeParams supplied directly. Useful for computing u_s for rotated shapes.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.scatteredfield-Tuple{Any,Any,ParticleScattering.ShapeParams,Any}",
    "page": "API",
    "title": "ParticleScattering.scatteredfield",
    "category": "method",
    "text": "scatteredfield(sigma_mu, k, s::ShapeParams, p) -> u_s\n\nComputes field scattered by the particle s with pre-computed potential densities sigma_mu at points p. All points must either be inside k = kin or outside k = kout the particle.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.solve_particle_scattering",
    "page": "API",
    "title": "ParticleScattering.solve_particle_scattering",
    "category": "function",
    "text": "solve_particle_scattering(k0, kin, P, sp::ScatteringProblem, θ_i = 0.0; get_inner = true, verbose = true) -> beta, inner\n\nSolve the scattering problem sp with outer wavenumber k0, inner wavenumber kin, 2P+1 cylindrical harmonics per inclusion and incident plane wave angle θ_i. Solves multiple-scattering equation directly. Returns the cylindrical harmonics basis beta along with potential densities (in case of arbitrary inclusion) or inner cylindrical coefficients (in case of circular). By default, incident wave propagates left->right.\n\nInner coefficients are only calculated if get_inner is true, and timing is printed if verbose is true.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.solve_particle_scattering_FMM-Tuple{Any,Any,Any,ParticleScattering.ScatteringProblem,Any,ParticleScattering.FMMoptions}",
    "page": "API",
    "title": "ParticleScattering.solve_particle_scattering_FMM",
    "category": "method",
    "text": "solve_particle_scattering_FMM(k0, kin, P, sp::ScatteringProblem, θ_i, opt::FMMoptions; plot_res = false, get_inner = true, verbose = true) -> result, inner\n\nSolve the scattering problem sp with outer wavenumber k0, inner wavenumber kin, 2P+1 cylindrical harmonics per inclusion and incident plane wave angle θ_i. Utilizes FMM with options opt to solve multiple-scattering equation. Returns the cylindrical harmonics basis beta along with convergence data in result. inner contains potential densities (in case of arbitrary inclusion) or inner cylindrical coefficients (in case of circular).\n\nplot_res controls plotting of the residual. Inner coefficients are calculated only if get_inner is true, and timing is printed if verbose is true.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.square_grid-Tuple{Integer,Any}",
    "page": "API",
    "title": "ParticleScattering.square_grid",
    "category": "method",
    "text": "square_grid(a::Integer, d)\n\nReturn centers, an (a^2,2) array containing the points on an a by a grid of points distanced d.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.squircle-Tuple{Any,Any}",
    "page": "API",
    "title": "ParticleScattering.squircle",
    "category": "method",
    "text": "squircle(r, N)\n\nReturn a ShapeParams object containing the shape parametrized by x()^4 + y()^4 = r^4 with 2N nodes.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.uniqueind-Union{Tuple{Array{T,1}}, Tuple{T}} where T<:Number",
    "page": "API",
    "title": "ParticleScattering.uniqueind",
    "category": "method",
    "text": "uniqueind(v::Vector{T}) where T <: Number -> inds,u\n\nGiven a vector of numbers v of length n, returns the unique subset u as well as a vector of indices inds of length n such that v == u[inds].\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.verify_min_distance-Tuple{Any,Array{Float64,2},Any,Array{Float64,2}}",
    "page": "API",
    "title": "ParticleScattering.verify_min_distance",
    "category": "method",
    "text": "verify_min_distance(shapes, centers::Array{Float64,2}, ids, points::Array{Float64,2})\nverify_min_distance(sp::ScatteringProblem, points)\n\nReturns true if the shapes placed at centers are properly distanced (non-intersecting scattering disks), and all points are outside the scattering disks.\n\n\n\n"
},

{
    "location": "api.html#ParticleScattering.verify_min_distance-Tuple{Any,Array{Float64,2},Any}",
    "page": "API",
    "title": "ParticleScattering.verify_min_distance",
    "category": "method",
    "text": "verify_min_distance(shapes, centers::Array{Float64,2}, ids)\nverify_min_distance(sp::ScatteringProblem)\n\nReturns true if the shapes placed at centers are properly distanced (non-intersecting scattering disks).\n\n\n\n"
},

{
    "location": "api.html#User-Interface-1",
    "page": "API",
    "title": "User Interface",
    "category": "section",
    "text": "Modules = [ParticleScattering]"
},

]}
