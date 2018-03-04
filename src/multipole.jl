"""
	solve_particle_scattering(k0, kin, P, sp::ScatteringProblem, θ_i = 0.0; get_inner = true, verbose = true) -> beta, inner

Solve the scattering problem `sp` with outer wavenumber `k0`, inner wavenumber
`kin`, `2P+1` cylindrical harmonics per inclusion and incident plane wave angle
`θ_i`. Solves multiple-scattering equation directly.
Returns the cylindrical harmonics basis `beta` along with potential
densities (in case of arbitrary inclusion) or inner cylindrical coefficients (in
case of circular). By default, incident wave propagates left->right.

Inner coefficients are only calculated if `get_inner` is true, and timing is
printed if `verbose` is true.
"""
function solve_particle_scattering(k0, kin, P, sp::ScatteringProblem, θ_i = 0.0;
								get_inner = true, verbose = true)
	# This function solves for the outgoing multipole coefficients in the presence
	# of an incident plane wave.
	# incident wave direction - from left to right is 0:
	# e^{ik(\cos(\θ_i),\sin(\θ_i)) \cdot \mathbf{r}}
	shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
	Ns = size(sp)
	#first solve for single scatterer densities
	α = Complex{Float64}[exp(1.0im*p*(pi/2-θ_i)) for p=-P:P]
	tic()
	scatteringMatrices,innerExpansions = particleExpansion(k0, kin, shapes, P, ids)
	dt1 = toq()
	tic()
	T = Array{Complex{Float64}}(Ns*(2*P+1),Ns*(2*P+1))
	a = Array{Complex{Float64}}(Ns*(2*P+1))
	Btemp = Array{Complex{Float64}}(2*P+1,2*P+1)
	for ic1 = 1:Ns
		rng1 = (ic1-1)*(2*P+1) + (1:2*P+1)
		if φs[ic1] == 0.0 || typeof(shapes[ids[ic1]]) == CircleParams
			RotScatMat = scatteringMatrices[ids[ic1]]
		else
			Rot = spdiagm(Complex{Float64}[exp(-1.0im*φs[ic1]*l) for l=-P:P]) #rotation matrix
			RotScatMat = Rot*(scatteringMatrices[ids[ic1]]*conj(Rot))
		end
		for ic2 = 1:Ns
			rng2 = (ic2-1)*(2*P+1) + (1:2*P+1)
			if ic1 == ic2
				T[rng1,rng2] = eye(Complex{Float64},2*P+1)
			else
				M2Lmatrix!(Btemp, k0, P, centers[ic1,:] - centers[ic2,:])
				T[rng1,rng2] = -RotScatMat*Btemp
			end
		end
		phase = exp(1.0im*k0*(cos(θ_i)*centers[ic1,1] + sin(θ_i)*centers[ic1,2])) #phase shift added to move cylinder coords
		a[rng1] = phase*(RotScatMat*α)
	end
	dt2 = toq()
	tic()
	beta = T\a
	dt3 = toq()
	#recover full incoming expansion - in sigma_mu terms for parametrized shape,
	#in multipole expansion for circle
	if get_inner
		tic()
		#find LU factorization once for each shape
        scatteringLU = [lufact(scatteringMatrices[iid]) for iid = 1:length(shapes)]
		inner = Array{Vector{Complex{Float64}}}(Ns)
		α_c = Array{Complex{Float64}}(2*P+1)
		for ic = 1:Ns
			rng = (ic-1)*(2*P+1) + (1:2*P+1)
			if typeof(shapes[ids[ic]]) == ShapeParams
				if φs[ic] == 0.0
					α_c[:] = scatteringLU[ids[ic]]\beta[rng]
				else
					Rot = spdiagm(Complex{Float64}[exp(-1.0im*φs[ic]*l) for l=-P:P]) #rotation matrix
					α_c[:] = scatteringLU[ids[ic]]\(conj(Rot)*beta[rng])
				end
				inner[ic] = innerExpansions[ids[ic]]*α_c
			else
				inner[ic] = innerExpansions[ids[ic]]*beta[rng]
			end
		end
		dt4 = toq()
	end
	verbose && begin
		println("Direct solution timing:")
		println("Scattering matrix solution: $dt1 s")
		println("Matrix construction: $dt2 s")
		println("Matrix solution: $dt3 s")
        get_inner && println("Retrieving inner coefficients: $dt4 s")
    end
	get_inner ? (return beta, inner) : (return beta)
end

function M2Lmatrix!(T, k, P, d)
	#builds non-diagonal translation matrix from particle 2 to particle 1 (d=1-2)
	#only calculates 2*(2*P+1) values
	kd = k*sqrt(sum(abs2,d))
	td = atan2(d[2],d[1])
	bess = besselh.(0:2*P,1,kd)
	for ix = 1:2*P #lower diagonals
		rng = ix+1:1+(2*P+1):(2*P+1)^2-(2*P+1)*ix
		T[rng] = exp(-1im*td*ix)*(-1)^(ix)*bess[ix+1]
	end
	for ix = 0:2*P #central and upper diagonals
		rng = ix*(2*P+1)+1:1+(2*P+1):(2*P+1)^2-ix
		T[rng] = exp(1im*td*ix)*bess[ix+1]
	end
end

function scattered_field_multipole(k0, beta, centers::Array{Float64,2}, points::Array{Float64,2})
	Ns = size(centers,1)
	P = div(div(length(beta),Ns)-1,2)
	len = size(points,1)
	Esc = zeros(Complex{Float64}, len)

	scattered_field_multipole!(Esc, k0, beta, P, centers, 1:Ns, points, 1:len)
	return Esc
end

function scattered_field_multipole!(Esc::Array{Complex{Float64},1}, k0, beta, P, centers::Array{Float64,2}, ind_centers, points::Array{Float64,2}, ind_points)
	for ic in ind_centers
		ind = (ic-1)*(2*P+1) + P + 1
		for ip in ind_points
			points_moved1 = points[ip,1] - centers[ic,1]
			points_moved2 = points[ip,2] - centers[ic,2]

			rs_moved = hypot(points_moved1, points_moved2)
			ts_moved = atan2(points_moved2, points_moved1)
			try
				Esc[ip] += beta[ind]*besselh(0, 1, k0*rs_moved)
			catch
				error("rs_moved == 0, center=$(centers[ic,:]), point=$(points[ip,:]), k0 = $k0, ic=$ic, ip=$ip")
			end
			for p = 1:P
				Esc[ip] += besselh(p, 1, k0*rs_moved)*(beta[p + ind]*exp(1im*p*ts_moved) + (-1)^p*beta[-p + ind]*exp(-1im*p*ts_moved))
			end
		end
	end
end

function circleScatteringMatrix(kout, kin, R, P; gamma = false)
    #non-vectorized, reuses bessel
    S = Array{Complex{Float64}}(2*P+1)
    gamma && (G = Array{Complex{Float64}}(2*P+1))

    pre_J0 = besselj(-P-1,kout*R)
    pre_J1 = besselj(-P-1,kin*R)
    pre_H = besselh(-P-1,kout*R)
    for p = -P:0
        J0 = besselj(p,kout*R)
        J1 = besselj(p,kin*R)
        H = besselh(p,kout*R)

        dJ0 = kout*(pre_J0 - (p/kout/R)*J0)
        dJ1 = kin*(pre_J1 - (p/kin/R)*J1)
        dH = kout*(pre_H - (p/kout/R)*H)

        S[p+P+1] = -(dJ0*J1 - J0*dJ1)/(dH*J1 - H*dJ1)
        p != 0 && (S[P+1-p] = S[p+P+1])

        if gamma #also returns gamma_n as a function of beta_n
            G[p+P+1] = (J0*dH - dJ0*H)/(J0*dJ1 - dJ0*J1)
            p != 0 && (G[P+1-p] = G[p+P+1])
        end
        pre_J0 = J0
        pre_J1 = J1
        pre_H = H
    end
    ScatMat = spdiagm(S)
    if gamma
        InnerMat = spdiagm(G)
        return ScatMat,InnerMat
    else
        return ScatMat
    end
end

function innerFieldCircle(kin, gamma, center::Array{Float64,1}, points::Array{Float64,2})
    #non-vectorized, reuses bessel
    P = div(length(gamma)-1,2)

	len = size(points,1)
	points_moved = similar(points)
    points_moved[:,1] = points[:,1] - center[1]
	points_moved[:,2] = points[:,2] - center[2]

	rs_moved = sqrt.(sum(abs2,points_moved,2))
	ts_moved = atan2.(points_moved[:,2], points_moved[:,1])

	bess = [besselj(p,kin*rs_moved[ii]) for ii=1:len, p=0:P]
	E = gamma[P + 1]*bess[:,1]
	for p = 1:P
		E += bess[:,p+1].*(gamma[p + P + 1]*exp.(1.0im*p*ts_moved) + (-1)^p*gamma[-p + P + 1]*exp.(-1.0im*p*ts_moved))
	end
	return E
end

function innerFieldCircle(kin, gamma, center::Array{Float64,1}, points::Array{Float64,1})
    #non-vectorized, reuses bessel
    P = div(length(gamma)-1,2)

    points_moved = points - center
	rs_moved = sqrt(sum(abs2,points_moved))
	ts_moved = atan2(points_moved[2], points_moved[1])

	bess = [besselj(p,kin*rs_moved) for p=0:P]
	E = gamma[P + 1]*bess[1]
	for p = 1:P
		E += bess[p+1]*(gamma[p + P + 1]*exp(1.0im*p*ts_moved) + (-1)^p*gamma[-p + P + 1]*exp(-1.0im*p*ts_moved))
	end
	return E
end

function particleExpansion(k0, kin, shapes, P, ids)
    scatteringMatrices = Array{Union{SparseMatrixCSC{Complex{Float64},Int64},Array{Complex{Float64},2}}}(0)
	innerExpansions = Array{Union{SparseMatrixCSC{Complex{Float64},Int64},Array{Complex{Float64},2}}}(0)
	for i = 1:length(shapes)
        #no use in computing matrices if shape doesn't actually show up! push garbage to maintain order
		if all(ids .!= i)
			push!(scatteringMatrices, [0.0im 0.0im])
			push!(innerExpansions, [0.0im 0.0im])
			continue
		end
		if typeof(shapes[i]) == ShapeParams
            AB = shapeMultipoleExpansion(k0, shapes[i].t, shapes[i].ft, shapes[i].dft, P)
    		sigma_mu_mult = get_potential(k0, kin, P, shapes[i].t, shapes[i].ft, shapes[i].dft)
    		ScatMat = AB*sigma_mu_mult
            push!(scatteringMatrices, ScatMat)
            push!(innerExpansions, sigma_mu_mult)
		else
			ScatMat, inner = circleScatteringMatrix(k0, kin, shapes[i].R, P, gamma = true)
			push!(scatteringMatrices, ScatMat)
			push!(innerExpansions, inner)
		end
	end
    return scatteringMatrices,innerExpansions
end

function rotateMultipole!(arr,phi,P)
    #applies e^{-i \phi l}, p=-P...P
    #if slow, change from view to direct indexing (should be negligible despite higher memory use)
    for p = -P:P
        arr[P + 1 + p] *= exp(-1.0im*phi*p)
    end
end

function rotateMultipole!(dest_arr,source_arr,phi,P)
    #applies e^{-i \phi l}, p=-P...P
    #if slow, change from view to direct indexing (should be negligible despite higher memory use)
    for p = -P:P
        dest_arr[P + 1 + p] = source_arr[P + 1 + p]*exp(-1.0im*phi*p)
    end
end
