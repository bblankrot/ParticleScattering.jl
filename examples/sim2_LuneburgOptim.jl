using PyPlot, ParticleScattering
import Optim, JLD

output_dir = "/home/bblankro/Dropbox/outputdir"
er = 4.5
k0 = 2π
l0 = 2π/k0
a_lens = 0.2*l0
R_lens = 10*a_lens

kin = k0*sqrt(er)
N_cells = Int64(round(2*R_lens/a_lens))
centers, ids_lnbrg, rs_lnbrg = luneburg_grid(R_lens, N_cells, er)
φs = zeros(Float64, length(ids_lnbrg))
θ_i = 0.0
P = 5

fmm_options = FMMoptions(true, acc = 6, dx = 2*a_lens, method = "pre")

optim_options =  Optim.Options(x_tol = 1e-6, iterations = 100, store_trace = true, extended_trace = true, show_trace = true, allow_f_increases = true)
linesearch = LineSearches.BackTracking()

points = [R_lens 0.0]
r_max = (a_lens/1.15/2)*ones(size(centers,1))
r_min = (a_lens*1e-3)*ones(size(centers,1))
rs0 = (0.25*a_lens)*ones(size(centers,1))

ids_max = collect(1:length(rs0))
test_max = optimize_radius(rs0, r_min, r_max, points, ids_max, P, θ_i, k0, kin, #precompile
                centers, fmm_options, optim_options, false, "BFGS", linesearch)
tic()
test_max = optimize_radius(rs0, r_min, r_max, points, ids_max, P, θ_i, k0, kin,
                centers, fmm_options, optim_options, false, "BFGS", linesearch)
optim_time = toq()
rs_max = test_max.minimizer

# plot near fields
filename1 = output_dir * "/opt_r_luneburg.tex"
filename2 = output_dir * "/opt_r_max.tex"
filename3 = output_dir * "/opt_r_0.tex"
border = (R_lens + a_lens)*[-1;1;-1;1]

sp1 = ScatteringProblem([CircleParams(rs_lnbrg[i]) for i in eachindex(rs_lnbrg)],
        ids_lnbrg, centers, φs)
Ez1 = plot_near_field(k0, kin, P, sp1, θ_i, x_points = 150, y_points = 150,
        opt = fmm_options, border = border)
plotNearField_pgf(filename1, k0, kin, P, sp1, θ_i; opt = fmm_options,
    x_points = 201, y_points = 201, border = border)

sp2 = ScatteringProblem([CircleParams(rs_max[i]) for i in eachindex(rs_max)],
        ids_max, centers, φs)
Ez2 = plot_near_field(k0, kin, P, sp2, θ_i, x_points = 150, y_points = 150,
            opt = fmm_options, border = border)
plotNearField_pgf(filename2, k0, kin, P, sp2, θ_i; opt = fmm_options,
    x_points = 201, y_points = 201, border = border)

sp3 = ScatteringProblem([CircleParams(rs0[i]) for i in eachindex(rs0)],
        collect(1:length(rs0)), centers, φs)
Ez3 = plot_near_field(k0, kin, P, sp3, θ_i, x_points = 150, y_points = 150,
        opt = fmm_options, border = border)
plotNearField_pgf(filename3, k0, kin, P, sp3, θ_i; opt = fmm_options,
    x_points = 201, y_points = 201, border = border)

#plot convergence
inner_iters = length(test_max.trace)
iters = [test_max.trace[i].iteration for i=1:inner_iters]
fobj = -[test_max.trace[i].value for i=1:inner_iters]
gobj = [test_max.trace[i].g_norm for i=1:inner_iters]
rng = iters .== 0

test_max_trace = test_max.trace
trace_of_r = [test_max.trace[i].metadata["x"] for i=1:inner_iters]
JLD.@save output_dir * "/luneburg_optim.jld"

# figure()
# plot(0:inner_iters-1, fobj)
# plot(0:inner_iters-1, gobj)
# plot((0:inner_iters-1)[rng], fobj[rng],"*")
# plot((0:inner_iters-1)[rng], gobj[rng],"*")

import PGFPlotsX; const pgf = PGFPlotsX
pgf.@pgf begin
    fobj_plot = pgf.Plot(pgf.Coordinates(0:inner_iters-1, fobj),
                {blue, thick, no_markers},
                label = "\$f_{\\mathrm{obj}}\$")
    gobj_plot = pgf.Plot(pgf.Coordinates(0:inner_iters-1, gobj),
                {red, thick, dashed, no_markers},
                label = "\$\\|\\mathbf{g}_{\\mathrm{obj}}\\|_{\\infty}\$")
    fobj_outer = pgf.Plot(pgf.Coordinates((0:inner_iters-1)[rng], fobj[rng]),
                {blue, only_marks, mark = "*", mark_options = {fill = "blue"}})
    gobj_outer = pgf.Plot(pgf.Coordinates((0:inner_iters-1)[rng], gobj[rng]),
                {red, only_marks, mark = "triangle*", mark_options = {fill = "red"}})
    ax = pgf.Axis([fobj_plot;gobj_plot;fobj_outer;gobj_outer],
        {
            width = "\\figurewidth",
            xlabel = "Iterations",
            legend_pos = "north east",
            legend_style = "font = \\footnotesize"
        })
end
pgf.save(output_dir * "/opt_r_conv.tex", ax ,include_preamble = false)

################ Testing with symmetry ######################
assert(length(ids_max)==size(centers,1))
centers_abs = centers[:,1] + 1im*abs.(centers[:,2])
ids_sym, centers_abs = ParticleScattering.my_uniqueind(centers_abs)
J = length(centers_abs)
r_max = (a_lens/1.15/2)*ones(J)
r_min = (a_lens*1e-3)*ones(J)
rs0 = (0.25*a_lens)*ones(J)

tic()
test_max_sym = optimize_radius(rs0, r_min, r_max, points, ids_sym, P, θ_i, k0, kin,
                centers, fmm_options, optim_options, false, "BFGS", linesearch)
sym_time = toq()
rs_sym = test_max_sym.minimizer
JLD.@save output_dir * "/luneburg_optim_sym.jld" test_max_sym sym_time

sp4 = ScatteringProblem([CircleParams(rs_sym[i]) for i in eachindex(rs_sym)],
        ids_sym, centers, φs)
Ez4 = plot_near_field(k0, kin, P, sp4, θ_i, x_points = 150, y_points = 150,
        opt = fmm_options, border = border)

u1 = calc_near_field(k0, kin, 7, sp1, points, θ_i; opt = fmm_options)
u2 = calc_near_field(k0, kin, 7, sp2, points, θ_i; opt = fmm_options)
u3 = calc_near_field(k0, kin, 7, sp3, points, θ_i; opt = fmm_options)
u4 = calc_near_field(k0, kin, 7, sp4, points, θ_i; opt = fmm_options)
abs.([u1[1];u2[1];u3[1];u4[1]])

#####################################
# selfconsistent err P calculation
Ez_4 = calc_near_field(k0, kin, 4, sp1, points, θ_i; opt = fmm_options)
Ez_5 = calc_near_field(k0, kin, 5, sp1, points, θ_i; opt = fmm_options)
Ez_6 = calc_near_field(k0, kin, 6, sp1, points, θ_i; opt = fmm_options)
Ez_7 = calc_near_field(k0, kin, 7, sp1, points, θ_i; opt = fmm_options)
