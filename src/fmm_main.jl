"""
	solve_particle_scattering_FMM(k0, kin, P, sp::ScatteringProblem, pw::PlaneWave, opt::FMMoptions; plot_res = false, get_inner = true, verbose = true) -> result, inner

Solve the scattering problem `sp` with outer wavenumber `k0`, inner wavenumber
`kin`, `2P+1` cylindrical harmonics per inclusion and incident plane wave angle
`pw.θi`. Utilizes FMM with options `opt` to solve multiple-scattering equation.
Returns the cylindrical harmonics basis `beta` along with convergence data in
`result`. `inner` contains potential densities (in case of arbitrary inclusion)
or inner cylindrical coefficients (in case of circular).

`plot_res` controls plotting of the residual. Inner coefficients are calculated
only if `get_inner` is true, and timing is printed if `verbose` is true.
"""
#TODO: return beta,inner,history
function solve_particle_scattering_FMM(k0, kin, P, sp::ScatteringProblem, pw::PlaneWave, opt::FMMoptions; plot_res = false, get_inner = true, verbose = true)
	assert(opt.FMM)
	shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    Ns = size(sp)
    groups, boxSize = divideSpace(centers, opt)
    P2, Q = FMMtruncation(opt.acc, boxSize, k0)
    verbose && println("FMM solution timing:")
    tic()
    mFMM = FMMbuildMatrices(k0, P, P2, Q, groups, centers, boxSize, tri=true)
    dt0 = toq()

    tic()
    scatteringMatrices,innerExpansions = particleExpansion(k0, kin, shapes, P, ids)
    dt1 = toq()

    tic()
    #construct rhs
    rhs = u2α(k0, pw, centers, P)
    for ic = 1:Ns
        rng = (ic-1)*(2*P+1) + (1:2*P+1)
        #see if there is a faster alternative
        if φs[ic] == 0.0
            rhs[rng] = scatteringMatrices[ids[ic]]*rhs[rng]
        else
            #rotate without matrix
            rotateMultipole!(view(rhs,rng),-φs[ic],P)
            rhs[rng] = scatteringMatrices[ids[ic]]*rhs[rng]
            rotateMultipole!(view(rhs,rng),φs[ic],P)
        end
    end
    dt2 = toq()

    pre_agg_buffer = zeros(Complex{Float64},Q,length(groups))
    trans_buffer = Array{Complex{Float64}}(Q)
    tic()
    if opt.method == "pre"
        MVP = LinearMap{eltype(rhs)}((output_, x_) -> FMM_mainMVP_pre!(output_,
									x_, scatteringMatrices, φs, ids, P, mFMM,
									pre_agg_buffer, trans_buffer), Ns*(2*P+1),
									Ns*(2*P+1), ismutating = true)
    elseif opt.method == "pre2"
        MVP = LinearMap{eltype(rhs)}((output_, x_) -> FMM_mainMVP_pre2!(output_,
		 							x_, scatteringMatrices, φs, ids, P, mFMM,
									pre_agg_buffer, trans_buffer), Ns*(2*P+1),
									Ns*(2*P+1), ismutating = true)
    end
    result = gmres(MVP, rhs, restart = Ns*(2*P+1), tol = opt.tol, log = true) #no restart, preconditioning
    dt3 = toq()

    if get_inner
        #recover full incoming expansion - in sigma_mu terms for parametrized shape,
        #in multipole expansion for circle
        tic()

        #find LU factorization once for each shape
        scatteringLU = [lufact(scatteringMatrices[i]) for i = 1:length(shapes)]

        inner = Array{Vector{Complex{Float64}}}(Ns)
        α_c = Array{Complex{Float64}}(2*P+1)
    	for ic = 1:Ns
    		rng = (ic-1)*(2*P+1) + (1:2*P+1)
            if typeof(shapes[ids[ic]]) == ShapeParams
                if φs[ic] == 0.0
                    α_c[:] = scatteringLU[ids[ic]]\result[1][rng]
                else
                    rotateMultipole!(α_c,view(result[1],rng),-φs[ic],P)
                    A_ldiv_B!(scatteringLU[ids[ic]],α_c)
                end
        		inner[ic] = innerExpansions[ids[ic]]*α_c
            else
                inner[ic] = innerExpansions[ids[ic]]*result[1][rng]
            end
    	end
        dt4 = toq()
    end

    if plot_res
        residual_ = MVP*result[1] - rhs
        println("last residual is empirically: $(norm(residual_)), but result[2].residuals[end] = $(result[2].residuals[end])")
        figure()
        semilogy(result[2].residuals.')
    end
    if verbose
        println("FMM matrix construction: $dt0 s")
        println("Scattering matrix solution: $dt1 s")
        println("RHS construction: $dt2 s")
        println("Iterative process: $dt3 s")
        get_inner && println("Retrieving inner coefficients: $dt4 s")
    end
    get_inner ? (return result, inner) : (return result)
end
