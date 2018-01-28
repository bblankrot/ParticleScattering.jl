function optimize_radius_dep(rs0, r_min, r_max, points, P, θ_i, k0, kin, centers, fmmopts, optimopts, minimize, method, linesearch)
    #assumes all are circles with their own IDs.
    #TODO: separate ids from rs to allow symmetry.
    #setup FMM reusal
    Ns = size(centers,1)
    (groups, boxSize) = divideSpace(centers, fmmopts)

    (P2, Q) = FMMtruncation(fmmopts.acc, boxSize, k0)
    mFMM = FMMbuildMatrices(k0, P, P2, Q, groups, centers, boxSize, tri=true)

    ids = collect(1:Ns) #all different sizes
    scatteringMatrices = Array{SparseMatrixCSC{Complex{Float64},Int64}}(0)
    #allocate derivative
    dS_S = Array{SparseMatrixCSC{Complex{Float64},Int64}}(0)
    for ic = 1:Ns
        push!(scatteringMatrices,speye(2*P+1))
        push!(dS_S,speye(2*P+1))
    end

    buf = FMMbuffer(Ns,P,Q,length(groups))

    #stuff that is done once
    points_moved = Array{Float64}(2)
    H = Array{Complex{Float64}}(Ns*(2*P+1),size(points,1))
    for i = 1:size(points,1), ic = 1:Ns
        points_moved[1] = points[i,1] - centers[ic,1]
        points_moved[2] = points[i,2] - centers[ic,2]
        r_angle = atan2(points_moved[2],points_moved[1])
        kr = k0*sqrt(points_moved[1]^2 + points_moved[2]^2)
        for l = -P:P
            H[(ic-1)*(2*P+1) + P + 1 + l,i] = besselh(l,kr)*exp(1.0im*l*r_angle)
        end
    end
    α_inc = Complex{Float64}[exp(1.0im*p*(pi/2-θ_i)) for p=-P:P]
    φs = zeros(Float64,Ns)

    # Allocate buffers
    #TODO: accept vecotr and scalar r_min,r_max
    shared_var = OptimBuffer(Ns,P,size(points,1))
    if length(rs0) == 1
        initial_rs = rs0*ones(Float64,Ns)
    else
        initial_rs = copy(rs0)
    end
    last_rs = similar(initial_rs)

    if minimize
        df = OnceDifferentiable(rs -> optimize_radius_f_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers,scatteringMatrices, dS_S, ids, mFMM, fmmopts, buf),
                    (grad_stor, rs) -> optimize_radius_g!_dep(grad_stor, rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, fmmopts, buf, minimize))
    else
        df = OnceDifferentiable(rs -> -optimize_radius_f_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers,scatteringMatrices, dS_S, ids, mFMM, fmmopts, buf),
                    (grad_stor, rs) -> optimize_radius_g!_dep(grad_stor, rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, fmmopts, buf, minimize))
    end

    outer_iterations = optimopts.iterations #?
    inner_iterations = optimopts.iterations #?

    if method == "LBFGS"
        optimize(df, initial_rs, r_min, r_max, Fminbox{LBFGS}();
            optimizer_o = optimopts, iterations = outer_iterations,
            linesearch = linesearch, x_tol = optimopts.x_tol,
            f_tol = optimopts.f_tol, g_tol = optimopts.g_tol)
    elseif method == "BFGS"
        optimize(df, initial_rs, r_min, r_max, Fminbox{BFGS}();
            optimizer_o = optimopts, iterations = outer_iterations,
            linesearch = linesearch, x_tol = optimopts.x_tol,
            f_tol = optimopts.f_tol, g_tol = optimopts.g_tol)
    end
end

function optimize_radius_common!_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)
    if (rs != last_rs)
        copy!(last_rs, rs)
        #do whatever common calculations and save to shared_var
        #construct rhs
        for ic = 1:Ns
            rng = (ic-1)*(2*P+1) + (1:2*P+1)
            try
                updateCircleScatteringDerivative!(scatteringMatrices[ids[ic]], dS_S[ids[ic]], k0, kin, rs[ic], P)
            catch
                warn("Could not calculate derivatives for ic=$ic,k0=$k0,kin=$kin,R=$(rs[ic])")
                error()
            end
            buf.rhs[rng] = scatteringMatrices[ids[ic]]*α_inc
            #phase shift added to move cylinder coords
            phase = exp(1.0im*k0*(cos(θ_i)*centers[ic,1] + sin(θ_i)*centers[ic,2]))
            buf.rhs[rng] *= phase
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
        shared_var.β,ch = gmres!(shared_var.β, MVP, buf.rhs, restart = Ns*(2*P+1) + 1, maxiter = Ns*(2*P+1), tol = opt.tol, log = true) #no restart, preconditioning
        if !ch.isconverged
            error("FMM process did not converge")
        end
        shared_var.f[:] = H.'*shared_var.β
        shared_var.f[:] += exp.(1.0im*k0*(cos(θ_i)*points[:,1] + sin(θ_i)*points[:,2])) #incident
    end
end

function optimize_radius_f_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)
    optimize_radius_common!_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)

    func = sum(abs2,shared_var.f)
end

function optimize_radius_g!_dep(grad_stor, rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf, minimize)
    optimize_radius_common!_dep(rs, last_rs, shared_var, φs, α_inc, H, points, P, θ_i, Ns, k0, kin, centers, scatteringMatrices, dS_S, ids, mFMM, opt, buf)

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
    shared_var.rhs_grad[:] = 0.0
    shared_var.∂β[:] = 0.0
    # tempn = Array{Complex{Float64}}(2*P+1)
    for n = 1:Ns
        #compute n-th gradient
        rng = (n-1)*(2*P+1) + (1:2*P+1)
        shared_var.rhs_grad[rng] = dS_S[ids[n]]*shared_var.β[rng]
        shared_var.∂β[:,n], ch = gmres!(shared_var.∂β[:,n], MVP, shared_var.rhs_grad, restart = Ns*(2*P+1) + 1, maxiter = Ns*(2*P+1), tol = opt.tol, log = true)
        #warn("using dbdn_tol = 10*opt.tol = $(10*opt.tol)")

        if ch.isconverged == false
            display(ch)
            display("rs:")
            display(rs)
            display("β:")
            display(shared_var.β)
            display("rhs_grad:")
            display(shared_var.rhs_grad)
            display("∂β:")
            display(shared_var.∂β[:,n])
            error("FMM process did not converge for partial derivative $n/$Ns. ")
        end
        #prepare for next one
        shared_var.rhs_grad[rng] = 0.0
    end

    grad_stor[:] = ifelse(minimize,2,-2)*real(shared_var.∂β.'*(H*conj(shared_var.f)))
end

function updateCircleScatteringDerivative!(S, dS_S, kout, kin, R, P)
    #non-vectorized, reuses bessel. Find quicker way to update

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

		numer = (-2.0im/pi/R)*(kin^2 - kout^2)*J1^2
        denom = (dH*J1 - H*dJ1)

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
