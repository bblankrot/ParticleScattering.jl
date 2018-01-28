# Here is a simulation of solution time as a function of M - number of shapes (or sqrt of number of shapes)
using ParticleScattering, LinearMaps, IterativeSolvers
using PyPlot, LaTeXStrings

#loop definitions
sqrtM_vec = collect(5:30); M_vec = sqrtM_vec.^2
trials = 3
simlen = length(M_vec)
res_vec = Array{Float64}(simlen,trials)
iter_vec = Array{Float64}(simlen,trials)
mvp_vec = Array{Float64}(simlen,trials)
setup_vec = Array{Float64}(simlen,trials)

#variables
k0 = 10.0
kin = 1.5*k0
l0 = 2π/k0
a1 = 0.3*l0; a2 = 0.1*l0; a3 = 5
θ_i = 0.0
tol = 1e-6
dist = 0.9l0

myshapefun(N) = rounded_star(a1, a2, a3, N)
N = 342, errN = 9.97e-7
# N,errN = minimumN(k0, kin, myshapefun, tol = tol, N_points = 20_000)
shapes = [myshapefun(N)]
P = 10, errP = 9.57e-7
# P,errP = minimumP(k0, kin, shapes[1], tol = tol, N_points = 20_000,
#                             P_min = 1, P_max = 120)

tic()
scatteringMatrices,innerExpansions = particleExpansion(k0, kin, shapes, P, [1])
α = Complex{Float64}[exp(1.0im*p*(pi/2-θ_i)) for p=-P:P]
dt0 = toq()
for i=1:simlen-1
    tic()
    scatteringMatrices,innerExpansions = particleExpansion(k0, kin, shapes, P, [1])
    α = Complex{Float64}[exp(1.0im*p*(pi/2-θ_i)) for p=-P:P]
    dt0 += toq()
end
dt0 /= simlen


for is = 1:simlen, it = 1:trials
    #compute shape variables
    begin #setup
        sqrtM = sqrtM_vec[is]
        M = sqrtM^2
        centers = square_grid(sqrtM, dist)
        φs = rand(M)
        ids = ones(Int64,M)
        opt = FMMoptions(true, acc = Int(-log10(tol)), nx = div(sqrtM,2), method="pre")
    end
    tic()
    (groups, boxSize) = divideSpace(centers, opt)
    (P2, Q) = FMMtruncation(opt.acc, boxSize, k0)
    mFMM = FMMbuildMatrices(k0, P, P2, Q, groups, centers, boxSize, tri=true)

    #construct rhs
    rhs = repeat(α,outer=[M])
    for ic = 1:M
        rng = (ic-1)*(2*P+1) + (1:2*P+1)
        if φs[ic] == 0.0
            rhs[rng] = scatteringMatrices[ids[ic]]*α
        else
            #rotate without matrix
            ParticleScattering.rotateMultipole!(view(rhs,rng),-φs[ic],P)
            rhs[rng] = scatteringMatrices[ids[ic]]*rhs[rng]
            ParticleScattering.rotateMultipole!(view(rhs,rng),φs[ic],P)
        end
        #phase shift added to move cylinder coords
        phase = exp(1.0im*k0*(cos(θ_i)*centers[ic,1] + sin(θ_i)*centers[ic,2]))
        rhs[rng] .*= phase
    end
    pre_agg_buffer = zeros(Complex{Float64},Q,length(groups))
    trans_buffer = Array{Complex{Float64}}(Q)

    MVP = LinearMap{eltype(rhs)}((output_, x_) -> ParticleScattering.FMM_mainMVP_pre!(output_, x_, scatteringMatrices,
        φs, ids, P, mFMM, pre_agg_buffer, trans_buffer), M*(2*P+1), M*(2*P+1), ismutating = true)
    x = zero(rhs)
    setup_vec[is,it] = toq()

    tic()
    x,ch = gmres!(x, MVP, rhs, restart = M*(2*P+1), tol = opt.tol, log = true, initially_zero = true) #no restart, preconditioning
    res_vec[is,it] = toq()
    tic()
    rhs[:] = MVP*x
    mvp_vec[is,it] = toq()
    iter_vec[is,it] = ch.iters
end

#average over all simulations
res_vec = vec(mean(res_vec,2))
iter_vec = vec(mean(iter_vec,2))
mvp_vec = vec(mean(mvp_vec,2))
setup_vec = vec(mean(setup_vec,2))

a_total,b_total = linreg(log10.(M_vec), log10.(res_vec))
res_ana = (10^a_total)*(M_vec.^b_total)
a_mvp,b_mvp = linreg(log10.(M_vec), log10.(mvp_vec))
mvp_ana = (10^a_mvp)*(M_vec.^b_mvp)
a_setup,b_setup = linreg(log10.(M_vec), log10.(setup_vec))
setup_ana = (10^a_setup)*(M_vec.^b_setup)
semilogy(M_vec, res_vec,"bo")
semilogy(M_vec, res_ana, "b-")
semilogy(M_vec, mvp_vec, "k+")
semilogy(M_vec, mvp_ana, "k-")
semilogy(M_vec, setup_vec, "r^")
semilogy(M_vec, setup_ana, "r-")
legend(("Elapsed time (Sol.)", @sprintf("\$%fM^{%.2f}\$", 10^a_total, b_total),
        "Elapsed time (MVP.)", @sprintf("\$%fM^{%.2f}\$", 10^a_mvp, b_mvp),
        "Elapsed time (Setup)", @sprintf("\$%fM^{%.2f}\$", 10^a_setup, b_setup)), loc = "best")
xlabel("Number of Scatterers")

### plot with pgfplots
import PGFPlotsX; const pgf = PGFPlotsX
pgf.@pgf begin
    ax = pgf.Axis({xlabel = "Number of scatterers",
            ylabel = "\$ \\mathrm{Run time} \\ [\\mathrm{s}]\$",
            xmode = "linear",
            ymode = "log",
            width = "\\figurewidth",
            legend_pos = "north west",
            legend_style = "font = \\footnotesize"})
    push!(ax, pgf.Plot(pgf.Coordinates(M_vec, res_vec),
            {blue, "only marks", mark = "*"};
            label = "Elapsed time (solution)"))
            temp2 = floor(log10(10^a_total))
            temp1 = 10^a_total/10^temp2
    push!(ax, pgf.Plot(pgf.Coordinates(M_vec, res_ana),
            {blue, thick, no_markers};
            label = @sprintf("\$%.1f \\cdot 10^{%d} \\cdot M^{%.2f}\$",
                temp1, temp2, b_total)))
    push!(ax, pgf.Plot(pgf.Coordinates(M_vec, mvp_vec),
            {red, only_marks, mark = "triangle*", mark_options = {fill = "red"}};
            label= "Elapsed time (MVP)"))
            temp2 = floor(log10(10^a_mvp))
            temp1 = 10^a_mvp/10^temp2
    push!(ax, pgf.Plot(pgf.Coordinates(M_vec, mvp_ana),
            {red, thick, dashed, no_markers};
            label = @sprintf("\$%.1f \\cdot 10^{%d} \\cdot M^{%.2f}\$",
                temp1, temp2, b_mvp)))
end

pgf.save(dirname(@__FILE__) * "/sim1.tex", ax ,include_preamble = false)
