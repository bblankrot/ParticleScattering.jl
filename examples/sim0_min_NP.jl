######################################
## comparing analytical (fictitious sources) and computed fields
using PyPlot,ParticleScattering
import JLD

function unique_subset(N,v,rev = true)
	all(sort(v, rev = rev) .== v) || error("Only handles sorted vectors.")
	v_u = [v[1]]
	N_u = [N[1]]
	prev = v[1]
	for i = 2:length(N)
		if v[i] !== prev
			#encountered new val
			push!(v_u, v[i])
			push!(N_u, N[i])
			prev = v[i]
		end
	end
    N_u,v_u
end

k0 = 10.0
kin = 1.5*k0
l0 = 2π/k0
a1 = 0.3l0
a2 = 0.1l0
N_points = 20_000
θ_i = 0.0

myshapefun5(N) = rounded_star(a1,a2,5,N)
myshapefun_squircle(N) = squircle(a1+0.5a2,N)

shape_functions = [myshapefun5; myshapefun_squircle]

if !isfile(dirname(@__FILE__) * "/mindata.jld")
    Nvec5 = [collect(20:10:90);collect(100:20:6000)]
    Nvecs = [collect(20:10:90);collect(100:20:6000)]
    errNvec5 = Array{Float64}(length(Nvec5))
    errNvecs = Array{Float64}(length(Nvecs))

    s = myshapefun5(400) #just for radius
    err_points = [s.R*f(i*2*pi/N_points) for i=0:(N_points-1), f in (cos,sin)]
    E_ana = (0.25im*besselh(0,1,k0*s.R))*ones(Complex{Float64},N_points)
    E_comp = Array{Complex{Float64}}(length(E_ana))

    for i in eachindex(Nvec5)
        errNvec5[i] = ParticleScattering.minimumN_helper(Nvec5[i], k0, kin,
						myshapefun5, err_points, E_comp, E_ana)
        display("Nvec5: done with $i/$(length(Nvec5))")
    end

    s = myshapefun_squircle(400) #just for radius
    err_points = [s.R*f(i*2*pi/N_points) for i=0:(N_points-1), f in (cos,sin)]
    E_ana = (0.25im*besselh(0,1,k0*s.R))*ones(Complex{Float64},N_points)
    E_comp = Array{Complex{Float64}}(length(E_ana))

    for i in eachindex(Nvecs)
        errNvecs[i] = ParticleScattering.minimumN_helper(Nvecs[i], k0, kin,
						myshapefun_squircle, err_points, E_comp, E_ana)
        display("Nvecs: done with $i/$(length(Nvecs))")
    end
    N5, errN5 = unique_subset(Nvec5, errNvec5)
    Ns, errNs = unique_subset(Nvecs, errNvecs)
    JLD.@save (dirname(@__FILE__) * "/mindata.jld") k0 kin l0 a1 a2 N5 errN5 Ns errNs
else
    JLD.@load (dirname(@__FILE__) * "/mindata.jld")
end


###############################3

inds = findfirst(errNs .< 5e-10)
ind5 = findfirst(errN5 .< 5e-10)
N5 = N5[1:ind5]
Ns = Ns[1:inds]
errN5 = errN5[1:ind5]
errNs = errNs[1:inds]

# here binary search is suboptimal since we know that new P is "close" to old one, frequently P or P+1
function findMinP(N, errN, shapefun, N_points, k0, kin; P_last = 1, P_max = 100)
	E_multipole = Array{Complex{Float64}}(N_points)
	errP = zeros(Float64, length(errN))
	Pmin = zeros(Int64, length(errN))
	for iN in eachindex(errN)
		s = shapefun(N[iN])

		err_points = [2.0*s.R*f(i*2*pi/N_points) for i=0:(N_points-1), f in (cos,sin)]
		#compute direct solution for comparison
		inner = solvePotentialShapePW(k0, kin, s.t, s.ft, s.dft, 0.0)
		E_quadrature = scatteredField(inner, k0, s.t, s.ft, s.dft, err_points)
		for P = P_last:P_max
			err = ParticleScattering.minimumP_helper(k0, kin, s, P, N_points,
					err_points, E_quadrature, E_multipole)
			if err <= errN[iN]
				errP[iN] = err
				Pmin[iN] = P
				P_last = P
				break
			elseif P == P_max
				warn("Failed to find P for iN = $iN")
				return Pmin, errP
			end
		end
		display("findMinP: done with iN=$iN/$(length(errN)), matched P=$(Pmin[iN]) with $(errP[iN]) ≦ $(errN[iN])")
	end
	Pmin, errP
end

Ps,errPs = findMinP(Ns, errNs, myshapefun_squircle, N_points, k0, kin)
JLD.@save dirname(@__FILE__) * "/mindata_s_new.jld" Ns errNs Ps errPs

P5,errP5 = findMinP(N5, errN5, myshapefun5, N_points, k0, kin)
JLD.@save dirname(@__FILE__) * "/mindata_5_new.jld" N5 errN5 P5 errP5


##################################################3
#plot with pgfplots
import PGFPlotsX; const pgf = PGFPlotsX

JLD.@load dirname(@__FILE__) * "/mindata_5_new.jld"
JLD.@load dirname(@__FILE__) * "/mindata_s_new.jld"

pgf.@pgf begin
    N_plot5 = pgf.Plot(pgf.Coordinates(errN5,N5),
        {blue, dashdotdotted, no_markers, thick},
		label = "\$N_{\\mathrm{min}} \\, \\mathrm{(star)}\$")
    P_plot5 = pgf.Plot(pgf.Coordinates(errN5,P5),
        {"green!50!black", no_markers, thick},
		label = "\$P_{\\mathrm{min}} \\, \\mathrm{(star)}\$")
    N_plot = pgf.Plot(pgf.Coordinates(errNs,Ns),
        {red, dotted, no_markers, thick},
		label = "\$N_{\\mathrm{min}} \\, \\mathrm{(squircle)}\$")
    P_plot = pgf.Plot(pgf.Coordinates(errNs,Ps),
        {black, dashed, no_markers, thick},
		label = "\$P_{\\mathrm{min}} \\, \\mathrm{(squircle)}\$")
    ax = pgf.Axis([N_plot5, N_plot, P_plot5, P_plot],
        {
            xmin = 5e-10,
            xmax = 1e-1,
			xlabel = "\$\\Delta u\$",
            xmode = "log",
            ymode = "log",
            grid = "both",
            legend_pos = "north east",
            legend_style = "font = \\footnotesize",
			legend_cell_align = "left"
        })
end
pgf.save(dirname(@__FILE__) * "/tikz/minNP.tex", ax ,include_preamble = false)

shape = myshapefun5(200)
pgf_shape = pgf.Coordinates([shape.ft[:,1];shape.ft[1,1]],
                            [shape.ft[:,2];shape.ft[1,2]])
pgf.@pgf begin
    ax2 = pgf.Axis(pgf.Plot(pgf_shape, {thick, black, no_markers, fill = "black!20"}),
        {   axis_equal,
            ticks = "none"})
    Rd = shape.R
    push!(ax2, "\\addplot [black, dashed, thick, domain=0:2*pi,samples=100]({$(shape.R)*cos(deg(x))},{$(shape.R)*sin(deg(x))});")
    push!(ax2, "\\node at (axis cs: -0.28,-0.15) {\$D\$};")
    push!(ax2, "\\node at (axis cs: 0.25,-0.25) {\$k_0\$};")
    push!(ax2, "\\node at (axis cs: 0.0,0.0) {\$k_1\$};")
    push!(ax2, "\\node at (axis cs: 0.15,0.1) {\$\\partial \\Omega_1\$};")
end
pgf.save(dirname(@__FILE__) * "/tikz/inclusion.tex", ax2 ,include_preamble = false)


shape = myshapefun_squircle(200)
pgf_shape = pgf.Coordinates([shape.ft[:,1];shape.ft[1,1]],
                            [shape.ft[:,2];shape.ft[1,2]])
pgf.@pgf begin
    ax3 = pgf.Axis(pgf.Plot(pgf_shape, {thick, black, no_markers, fill = "black!20"}),
        {   axis_equal,
            ticks = "none"})
    Rd = shape.R
    push!(ax3, "\\addplot [black, dashed, thick, domain=0:2*pi,samples=100]({$(shape.R)*cos(deg(x))},{$(shape.R)*sin(deg(x))});")
    push!(ax3, "\\node at (axis cs: -0.28,-0.15) {\$D\$};")
    push!(ax3, "\\node at (axis cs: 0.25,-0.25) {\$k_0\$};")
    push!(ax3, "\\node at (axis cs: 0.0,0.0) {\$k_1\$};")
    push!(ax3, "\\node at (axis cs: 0.15,0.1) {\$\\partial \\Omega_1\$};")
end
pgf.save(dirname(@__FILE__) * "/tikz/inclusion2.tex", ax3 ,include_preamble = false)
