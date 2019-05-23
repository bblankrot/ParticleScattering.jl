using PyPlot, ParticleScattering
import JLD, Optim, PGFPlotsX; const pgf = PGFPlotsX

output_dir = homedir()
Ns = 100
k0 = 2π
kin = 3*k0
l0 = 2π/k0
a1=0.3*l0; a2=0.1*l0; a3=5;
dmin = R_multipole*2*(a1+a2)
θ_i = 0.5π

size_factor = 7
width = size_factor*3*l0
height = size_factor*l0
myshapefun(N) = rounded_star(a1,a2,a3,N)
points = [range(0.0, stop=width, length=10) height*ones(10)]

if !isfile(dirname(@__FILE__) * "/sim3data.jld")
    centers = randpoints(Ns, dmin, width, height, points)
else
    import JLD; JLD.@load(dirname(@__FILE__) * "/sim3data.jld", centers)
end

N,errN = (934, 9.97040926753751e-7) #
#N,errN = minimumN(k0,kin,myshapefun, tol = 1e-6, N_points = 20_000)
shapes = [myshapefun(N)]
P,errP = (12, 8.538711552646218e-7)#
#P,errP = minimumP(k0, kin, shapes[1], tol = 1e-6, N_points = 20_000, P_min = 1, P_max = 120)

φs = zeros(Float64,Ns)
ids = ones(Int, Ns)

fmm_options = FMMoptions(true, acc = 6, nx = 9, method="pre")

divideSpace(centers, fmm_options; drawGroups = false)

draw_fig = false

# verify and draw
begin
    assert(verify_min_distance(shapes, centers, ids, points))
    if draw_fig
        figure()
        #draw shapes and points
        draw_shapes(shapes, ids, centers, φs)
        plot(points[:,1], points[:,2], "r*")
        tight_layout()
        ax = gca()
        ax[:set_aspect]("equal", adjustable = "box")
    end
end

plot_border = shapes[1].R*[-1;1;-1;1] + [0.0; width; 0.0; height]

optim_options =  Optim.Options(f_tol = 1e-6,
                                iterations = 100,
                                store_trace = true,
                                extended_trace = false,
                                show_trace = true,
                                allow_f_increases = true)

optim_method = Optim.BFGS(;linesearch = LineSearches.BackTracking())

tic()
test_max = optimize_φ(φs, points, P, PlaneWave(θ_i), k0, kin, shapes,
            centers, ids, fmm_options, optim_options, optim_method; minimize = false)
optim_time = toq()

# %%

sp_before = ScatteringProblem(shapes, ids, centers, φs)
plot_near_field(k0, kin, P, sp_before, PlaneWave(θ_i),
                x_points = 600, y_points = 200, border = plot_border);
colorbar()
clim([0;5])
plot_near_field_pgf(output_dir * "/opt_phi_before.tex", k0, kin, P,
    sp_before, PlaneWave(θ_i); opt = fmm_options, x_points = 150, y_points = 50,
    border = plot_border, downsample = 10, include_preamble = true)

sp_after = ScatteringProblem(shapes, ids, centers, test_max.minimizer)
plot_near_field(k0, kin, P, sp_after, PlaneWave(θ_i),
                x_points = 600, y_points = 200, border = plot_border)
colorbar()
clim([0;5])
plot_near_field_pgf(output_dir * "/opt_phi_after.tex", k0, kin, P,
    sp_after, PlaneWave(θ_i); opt = fmm_options, x_points = 150, y_points = 50,
    border = plot_border, downsample = 10, include_preamble = true)

inner_iters = length(test_max.trace)
fobj = -[test_max.trace[i].value for i=1:inner_iters]
gobj = [test_max.trace[i].g_norm for i=1:inner_iters]

figure()
plot(0:inner_iters-1, fobj)
plot(0:inner_iters-1, gobj)

pgf.@pgf begin
    fobj_plot = pgf.Plot({blue, thick, no_markers},
                        pgf.Coordinates(0:inner_iters-1, fobj))
    gobj_plot = pgf.Plot({red, thick, dashed, no_markers},
                        pgf.Coordinates(0:inner_iters-1, gobj))
    ax = pgf.Axis({ width = "6cm",
                    xlabel = "Iterations",
                    legend_pos = "north east",
                    legend_style = "font = \\footnotesize"},
                fobj_plot, gobj_plot)
    push!(ax, pgf.Legend(["\$f_{\\mathrm{obj}}\$";
                        "\$\\|\\mathbf{g}_{\\mathrm{obj}}\\|\$"]))
end
pgf.save(output_dir * "/opt_phi_conv.tex", ax ,include_preamble = false)
