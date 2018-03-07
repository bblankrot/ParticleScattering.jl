"""
    optimize_radius(rs0, r_min, r_max, points, ids, P, θ_i, k0, kin, centers,
        fmmopts, optimopts::Optim.Options; minimize = true, method = "BFGS")

Optimize the radii of circular particles for minimization or maximization of the
field intensity at `points`, depending on `minimize`. Uses `Optim`'s `Fminbox`
box-contrained optimization to contain radii in feasible rangle, given in scalar
or vector form by `r_min` and `r_max`.

Here, `ids` allows for grouping particles - for example, to maintain symmetry of
the optimized device.
`optimopts` defines the convergence criteria and other optimization parameters
for both the inner and outer iterations. `method` can be either `"BFGS"` or
`"LBFGS"`. See the `Optim.Fminbox` documentation for more details.

Returns an object of type `Optim.MultivariateOptimizationResults`.
"""
function optimize_radius(rs0, r_min, r_max, points, ids, P, θ_i, k0, kin,
                        centers, fmmopts, optimopts::Optim.Options;
                        minimize = true, method = "BFGS")
    Ns = size(centers,1)
    J = length(rs0)

    assert(maximum(ids) <= J)
    if length(r_min) == 1
        r_min = r_min*ones(Float64,J)
    else
        assert(J == length(r_min))
    end
    if length(r_max) == 1
        r_max = r_max*ones(Float64,J)
    else
        assert(J == length(r_max))
    end
    verify_min_distance([CircleParams(r_max[i]) for i = 1:J], centers, ids,
        points) || error("Particles are too close.")

    #setup FMM reusal
    groups, boxSize = divideSpace(centers, fmmopts)
    P2, Q = FMMtruncation(fmmopts.acc, boxSize, k0)
    mFMM = FMMbuildMatrices(k0, P, P2, Q, groups, centers, boxSize, tri = true)

    #allocate derivative
    scatteringMatrices = [speye(Complex{Float64}, 2*P+1) for ic = 1:J]
    dS_S = [speye(Complex{Float64}, 2*P+1) for ic = 1:J]

    #stuff that is done once
    H = optimizationHmatrix(points, centers, Ns, P, k0)
    α_inc = Complex{Float64}[exp(1.0im*p*(π/2-θ_i)) for p=-P:P]
    φs = zeros(Float64,Ns)

    # Allocate buffers
    buf = FMMbuffer(Ns,P,Q,length(groups))
    shared_var = OptimBuffer(Ns,P,size(points,1),J)
    initial_rs = copy(rs0)
    last_rs = similar(initial_rs)

    if minimize
        df = OnceDifferentiable(rs -> optimize_radius_f(rs, last_rs, shared_var,
                                        φs, α_inc, H, points, P, θ_i, Ns, k0,
                                        kin, centers,scatteringMatrices, dS_S,
                                        ids, mFMM, fmmopts, buf),
                    (grad_stor, rs) -> optimize_radius_g!(grad_stor, rs, last_rs,
                                        shared_var, φs, α_inc, H, points, P, θ_i,
                                        Ns, k0, kin, centers, scatteringMatrices,
                                        dS_S, ids, mFMM, fmmopts, buf, minimize),
                                        initial_rs)
    else
        df = OnceDifferentiable(rs -> -optimize_radius_f(rs, last_rs, shared_var,
                                        φs, α_inc, H, points, P, θ_i, Ns, k0,
                                        kin, centers,scatteringMatrices, dS_S,
                                        ids, mFMM, fmmopts, buf),
                    (grad_stor, rs) -> optimize_radius_g!(grad_stor, rs, last_rs,
                                        shared_var, φs, α_inc, H, points, P, θ_i,
                                        Ns, k0, kin, centers, scatteringMatrices,
                                        dS_S, ids, mFMM, fmmopts, buf, minimize),
                                        initial_rs)
    end

    outer_iterations = optimopts.iterations

    if method == "LBFGS"
        optimize(df, initial_rs, r_min, r_max, Fminbox{LBFGS}();
            optimizer_o = optimopts, iterations = outer_iterations,
            linesearch = LineSearches.BackTracking(), x_tol = optimopts.x_tol,
            f_tol = optimopts.f_tol, g_tol = optimopts.g_tol)
    elseif method == "BFGS"
        optimize(df, initial_rs, r_min, r_max, Fminbox{BFGS}();
            optimizer_o = optimopts, iterations = outer_iterations,
            linesearch = LineSearches.BackTracking(), x_tol = optimopts.x_tol,
            f_tol = optimopts.f_tol, g_tol = optimopts.g_tol)
    end
end

function optimize_radius_common!(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)
    if (rs != last_rs)
        copy!(last_rs, rs)
        #do whatever common calculations and save to shared_var
        #construct rhs
        for id in ids
            try
                updateCircleScatteringDerivative!(scatteringMatrices[id], dS_S[id], k0, kin, rs[id], P)
            catch
                warn("Could not calculate derivatives for id=$id,k0=$k0,kin=$kin,R=$(rs[id])")
                error()
            end
        end
        for ic = 1:Ns
            rng = (ic-1)*(2*P+1) + (1:2*P+1)
            buf.rhs[rng] = scatteringMatrices[ids[ic]]*α_inc
            #phase shift added to move cylinder coords
            phase = exp(1.0im*k0*(cos(θ_i)*centers[ic,1] + sin(θ_i)*centers[ic,2]))
            buf.rhs[rng] .*= phase
        end

        if opt.method == "pre"
            MVP = LinearMap{eltype(buf.rhs)}((output_, x_) -> FMM_mainMVP_pre!(output_,
                x_, scatteringMatrices, φs, ids, P, mFMM, buf.pre_agg, buf.trans),
                Ns*(2*P+1), Ns*(2*P+1), ismutating = true)
        elseif opt.method == "pre2"
            MVP = LinearMap{eltype(buf.rhs)}((output_, x_) -> FMM_mainMVP_pre2!(output_,
                x_, scatteringMatrices, φs, ids, P, mFMM, buf.pre_agg, buf.trans),
                Ns*(2*P+1), Ns*(2*P+1), ismutating = true)
        end

        fill!(shared_var.β,0.0)
        #no restart
        shared_var.β,ch = gmres!(shared_var.β, MVP, buf.rhs,
                            restart = Ns*(2*P+1) + 1, maxiter = Ns*(2*P+1),
                            tol = opt.tol, log = true, initially_zero = true)
        if !ch.isconverged
            error("FMM process did not converge")
        end
        shared_var.f[:] = H.'*shared_var.β
        shared_var.f[:] += exp.(1.0im*k0*(cos(θ_i)*points[:,1] + sin(θ_i)*points[:,2])) #incident
    end
end

function optimize_radius_f(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)
    optimize_radius_common!(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)

    func = sum(abs2,shared_var.f)
end

function optimize_radius_g!(grad_stor, rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf, minimize)
    optimize_radius_common!(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)

    if opt.method == "pre"
        MVP = LinearMap{eltype(buf.rhs)}((output_, x_) -> FMM_mainMVP_pre!(output_,
            x_, scatteringMatrices, φs, ids, P, mFMM, buf.pre_agg, buf.trans),
            Ns*(2*P+1), Ns*(2*P+1), ismutating = true)
    elseif opt.method == "pre2"
        MVP = LinearMap{eltype(buf.rhs)}((output_, x_) -> FMM_mainMVP_pre2!(output_,
            x_, scatteringMatrices, φs, ids, P, mFMM, buf.pre_agg, buf.trans),
            Ns*(2*P+1), Ns*(2*P+1), ismutating = true)
    end

    #time for gradient
    shared_var.∂β[:] = 0.0
    for n = 1:length(rs)
        #compute n-th gradient - here we must pay the price for symmetry
        #as more than one beta is affected.
        #TODO:expand to phi optimization?, and clean+speed this up.
        for ic = 1:Ns
            rng = (ic-1)*(2*P+1) + (1:2*P+1)
            if ids[ic] == n
                shared_var.rhs_grad[rng] = dS_S[n]*shared_var.β[rng]
            else
                shared_var.rhs_grad[rng] = 0.0
            end
        end
        shared_var.∂β[:,n], ch = gmres!(shared_var.∂β[:,n], MVP,
                                    shared_var.rhs_grad,
                                    restart = Ns*(2*P+1) + 1,
                                    maxiter = Ns*(2*P+1), tol = opt.tol,
                                    log = true, initially_zero = true)
        if !ch.isconverged
            display(ch)
            display("rs:"); display(rs)
            display("β:"); display(shared_var.β)
            display("rhs_grad:"); display(shared_var.rhs_grad)
            display("∂β:"); display(shared_var.∂β[:,n])
            error("FMM process did not converge for partial derivative $n/$Ns.")
        end
    end

    grad_stor[:] = ifelse(minimize,2,-2)*real(shared_var.∂β.'*(H*conj(shared_var.f)))
end

function updateCircleScatteringDerivative!(S, dS_S, kout, kin, R, P)
    #non-vectorized, reuses bessel

    pre_J0 = besselj(-1,kout*R)
    pre_J1 = besselj(-1,kin*R)
    pre_H = besselh(-1,kout*R)
    for p = 0:P
        J0 = besselj(p,kout*R)
        J1 = besselj(p,kin*R)
        H = besselh(p,kout*R)

        dJ0 = kout*(pre_J0 - (p/kout/R)*J0)
        dJ1 = kin*(pre_J1 - (p/kin/R)*J1)
        dH = kout*(pre_H - (p/kout/R)*H)

		numer = (-2.0im/(π*R))*(kin^2 - kout^2)*J1^2
        denom = dH*J1 - H*dJ1

		S[p+P+1,p+P+1] = -(dJ0*J1 - J0*dJ1)/denom
		dS_S[p+P+1,p+P+1] = -(numer/denom)/(dJ0*J1 - J0*dJ1)

		if p != 0
			S[P+1-p,P+1-p] = S[p+P+1,p+P+1]
			dS_S[P+1-p,P+1-p] = dS_S[p+P+1,p+P+1]
		end

        pre_J0 = J0
        pre_J1 = J1
        pre_H = H
    end
end
