function SDNTpotentialsdiff(k1, k2, t, ft, dft)
    #now just returns system matrix, utilizes similarities between upper and lower triangles
    iseven(length(t)) ?
    (N = div(length(t),2)) : (error("length(t) must be even"))

    (Rvec, kmlogvec) = KM_weights(N)
    A = Array{Complex{Float64}}(4*N, 4*N)

    ndft = sqrt.(vec(sum(abs2,dft,2)))
    rij = Array{Float64}(2)
    for i=1:2*N, j=1:i
        if i == j
            T1 = -(k1^2 - k2^2)
            T2 = k1^2*(π*1im - 2γ + 1 - 2*log(k1*ndft[i]/2)) - k2^2*(π*1im - 2γ + 1 - 2*log(k2*ndft[i]/2))

            A[i,j] = (-ndft[i]/2/N)*log(k1/k2) #dS (dM1=0)
            A[i,j+2*N] = 1 #dD (dL=0)
            A[i+2*N,j] = -1 #dN=0
            A[i+2*N,j+2*N] = (ndft[i]/8π)*(Rvec[1]*T1 + (π/N)*T2) #dT
            continue
        end
        rij[1] = ft[i,1]-ft[j,1]
        rij[2] = ft[i,2]-ft[j,2] #ridiculous but much faster than ft[i,:]-ft[j,:]
        r = sqrt(rij[1]^2 +  rij[2]^2)
        didj = dft[i,1]*dft[j,1] + dft[i,2]*dft[j,2]

        J01 = besselj(0, k1*r)
        J02 = besselj(0, k2*r)
        H01 = besselh(0, 1, k1*r)
        H02 = besselh(0, 1, k2*r)
        k2J0 = k1^2*J01 - k2^2*J02
        k1J1 = k1*besselj(1,k1*r) - k2*besselj(1,k2*r)
        k2H0 = k1^2*H01 - k2^2*H02
        k1H1 = k1*besselh(1,k1*r) - k2*besselh(1,k2*r)

        kmlog = kmlogvec[abs(i-j)+1]
        R = Rvec[abs(i-j)+1]

        N2 = 1im*pi*k1H1 - k1J1*kmlog

        P1 = (-didj/π)*k2J0
        P2 = 1im*didj*k2H0 - P1*kmlog

        Qtilde = (dft[i,1]*rij[1] + dft[i,2]*rij[2])*(dft[j,1]*rij[1] + dft[j,2]*rij[2])/r^2

        Q1 = (-Qtilde*k2J0 + (1/r)*k1J1*(2*Qtilde - didj))/π
        Q2 = 1im*Qtilde*k2H0 + (-1im/r)*k1H1*(2*Qtilde - didj) - Q1*kmlog

        M1 = (J02 - J01)/π
        L1 = k1J1/π
        M2 = 1im*(H01 - H02) - M1*kmlog
        L2 = 1im*k1H1 - L1*kmlog

        wro_ij =  (dft[j,2]*rij[1] - dft[j,1]*rij[2])/r
        wro_ji = -(dft[i,2]*rij[1] - dft[i,1]*rij[2])/r
        cross_ij = -wro_ji*ndft[j]/ndft[i]
        cross_ji = -wro_ij*ndft[i]/ndft[j]

        #edited to remove division by wro which might be 0
        A[i,j] = (0.25*ndft[j])*(R*M1 + (π/N)*M2) #dS
        A[j,i] = A[i,j]*(ndft[i]/ndft[j]) #dS
        A[i,j+2*N] = (0.25*wro_ij)*(R*L1 + (π/N)*L2) #dD
        A[j,i+2*N] = (0.25*wro_ji)*(R*L1 + (π/N)*L2) #dD
        A[i+2*N,j] = (-0.25*cross_ij/π)*(R*k1J1 + (π/N)*N2) #dN
        A[j+2*N,i] = (-0.25*cross_ji/π)*(R*k1J1 + (π/N)*N2) #dN
        A[i+2*N,j+2*N] = (R*(P1-Q1) + (π/N)*(P2-Q2))/(4*ndft[i]) #dT
        A[j+2*N,i+2*N] = A[i+2*N,j+2*N]*(ndft[i]/ndft[j]) #dT
    end
    any(isnan.(A)) && error("SDNTpotentialsdiff: encountered NaN, check data and division by ndft.")
    return A
end

function solvePotentialShape(kout, kin, P, t, ft, dft)
    N = length(t) #N here is different...

    A = SDNTpotentialsdiff(kout, kin, t, ft, dft)
    LU = lufact(A) #w/o pivoting: (L,U) = lu(A,Val{false})

    sigma_mu = Array{Complex{Float64}}(2*N, 2*P+1)

    #assuming the wave is sampled on the shape
    nz = sqrt.(sum(abs2,ft,2))
    ndz = sqrt.(sum(abs2,dft,2))
    nzndz = nz.*ndz
    #precompute? or negative = (-1)^n positive?
    #no precompute
	wro = dft[:,2].*ft[:,1] - dft[:,1].*ft[:,2]
	zz = dft[:,1].*ft[:,1] + dft[:,2].*ft[:,2]

    bessp = besselj.(-P-1, kout*nz)
    bess = similar(bessp)
    du = Array{Complex{Float64}}(length(bessp))
    rhs = Array{Complex{Float64}}(2*length(bessp))
    for p = -P:P
        bess[:] = besselj.(p, kout*nz)
		du[:] = kout*bessp.*wro - (p*bess./nz).*(wro + 1im*zz)
        rhs[:] = -[bess.*exp.(1.0im*p*t);
               (du./nzndz).*exp.(1.0im*p*t)]
        sigma_mu[:,p + P + 1] = LU\rhs #w/o pivoting: sigma_mu[:,p + P + 1] = U\(L\rhs)
        copy!(bessp,bess)
    end
    return sigma_mu
end

function KM_weights(N)
    #computes the weights necessary for Kussmaul-Martensen quadrature (evenly
    #spaced).
    #Input: N (integer>=1)
    #Output: R,K (float vectors of length 2N)
    arg1 = Float64[cos(m*j*π/N)/m for m=1:N-1, j=0:2*N-1]
    R = vec((-2π/N)*sum(arg1,1)) - (π/N^2)*Float64[cos(j*π) for j=0:2*N-1]

    K = Float64[2*log(2*sin(0.5π*j/N)) for j = 0:2*N-1]
    return (R,K)
end

function solvePotentialShapePW(kout, kin, s, θ_i)
    N = length(s.t) #N here is different...

    A = SDNTpotentialsdiff(kout, kin, s.t, s.ft, s.dft)
    LU = lufact(A)

	ndft = sqrt.(sum(abs2,s.dft,2))
    uinc = exp.(1.0im*kout*(cos(θ_i)*s.ft[:,1] + sin(θ_i)*s.ft[:,2]))
    rhs = -[uinc;
			(1.0im*kout*uinc).*((cos(θ_i)*s.dft[:,2] - sin(θ_i)*s.dft[:,1])./ndft)]

    sigma_mu = LU\rhs
end

function scatteredField(sigma_mu, k, t, ft, dft, p)
    #calculates the scattered field of a shape with parametrization ft(t),...,dft(t)
    #in space with wavenumber k at points p *off* the boundary. For field on the boundary,
    #SDpotentials function must be used.
    N = length(t)
    M = size(p,1)
    #loop is faster here:
    SDout = Array{Complex{Float64}}(M, 2*N)
    for j = 1:N
        ndft = hypot(dft[j,1],dft[j,2])
        for i = 1:M
            r = [p[i,1] - ft[j,1];p[i,2] - ft[j,2]]
            nr = hypot(r[1],r[2])
            if nr < eps()
                #TODO: use SDNTpotentialsdiff here
                warn("Encountered singularity in scatteredField.")
                SDout[i,j] = 0
                SDout[i,j+N] = 0
                continue
            end
            SDout[i,j] = (2*pi/N)*0.25im*besselh(0,1, k*nr)*ndft
            SDout[i,j+N] = (2*pi/N)*0.25im*k*besselh(1,1, k*nr)*(dft[j,2]*r[1] - dft[j,1]*r[2])/nr
        end
    end
    u_scat = SDout*sigma_mu
end

function scatteredField(sigma_mu, kout, s::ShapeParams, p::Array{Float64,2})
    #calculates the scattered field of a shape with parametrization ft(t),...,dft(t)
    #in space with wavenumber kout at points p *off* the boundary. For field on the boundary,
    #SDpotentials function must be used.
    N = length(s.t)
    M = size(p,1)
    #loop is faster here:
    SDout = Array{Complex{Float64}}(M, 2*N)
    for j = 1:N
        ndft = hypot(s.dft[j,1], s.dft[j,2])
        for i = 1:M
            r = [p[i,1] - s.ft[j,1];p[i,2] - s.ft[j,2]]
            nr = hypot(r[1], r[2])
            if nr < eps()
                #TODO: use SDNTpotentialsdiff here
                warn("Encountered singularity in scatteredField.")
                SDout[i,j] = 0
                SDout[i,j+N] = 0
                continue
            end
            SDout[i,j] = besselh(0,1,kout*nr)*ndft
            SDout[i,j+N] = kout*besselh(1,1,kout*nr)*(s.dft[j,2]*r[1] - s.dft[j,1]*r[2])/nr
        end
    end
    u_scat = (0.5im*π/N)*(SDout*sigma_mu)
end

function scatteredField(sigma_mu, kout, s::ShapeParams, p::Array{Float64,1})
    #calculates the scattered field of a shape with parametrization ft(t),...,dft(t)
    #in space with wavenumber kout at points p *off* the boundary. For field on the boundary,
    #SDpotentials function must be used.
    N = length(s.t)
    SDout = Array{Complex{Float64}}(2*N)
    for j = 1:N
        ndft = hypot(s.dft[j,1], s.dft[j,2])
        i = 1
        r = [p[1] - s.ft[j,1]; p[2] - s.ft[j,2]]
        nr = hypot(r[1], r[2])
        if nr < eps()
            #TODO: use SDNTpotentialsdiff here
            warn("Encountered singularity in scatteredField.")
            SDout[j] = 0
            SDout[j+N] = 0
            continue
        end
        SDout[j] = besselh(0,1,kout*nr)*ndft
        SDout[j+N] = kout*besselh(1,1,kout*nr)*(s.dft[j,2]*r[1] - s.dft[j,1]*r[2])/nr
    end
    u_scat = (0.5im*π/N)*(SDout.'*sigma_mu)
end

function shapeMultipoleExpansion(k, t, ft, dft, P)
    #unlike others (so far), this does *not* assume t_j=pi*j/N
    N = div(length(t),2)
    nz = vec(sqrt.(sum(abs2,ft,2)))
    ndz = vec(sqrt.(sum(abs2,dft,2)))
    AB = Array{Complex{Float64}}(2*P + 1, 4*N)
    bessp = besselj.(-P-1,k*nz)
    bess = similar(bessp)
    for l = -P:0
        bess[:] = besselj.(l,k*nz)
        for j = 1:2*N
            AB[l+P+1,j] = 0.25im*(π/N)*bess[j]*exp(-1.0im*l*t[j])*ndz[j]
            l!=0 && (AB[-l+P+1,j] = 0.25im*((-1.0)^l*π/N)*bess[j]*exp(1.0im*l*t[j])*ndz[j])
            wro = ft[j,1]*dft[j,2] - ft[j,2]*dft[j,1]
            zdz = -1.0im*(ft[j,1]*dft[j,1] + ft[j,2]*dft[j,2])
            b1 = (-l*bess[j]/nz[j])*(zdz + wro)
            b1_ = (-l*bess[j]/nz[j])*(zdz - wro)
            b2 = k*bessp[j]*wro
            AB[l+P+1,j+2*N] = 0.25im*(π/N)*(exp(-1.0im*l*t[j])/nz[j])*(b1 + b2)
            l!=0 && (AB[-l+P+1,j+2*N] = 0.25im*((-1.0)^l*π/N)*(exp(1.0im*l*t[j])/nz[j])*(-b1_ + b2))
        end
        copy!(bessp,bess)
    end
    return AB
end

function solvePotential_forError(kin, kout, shape, ls_pos, ls_amp, θ_i)
    #plane wave outside, line sources inside
    N = length(shape.t) #N here is different...

    A = SDNTpotentialsdiff(kout, kin, shape.t, shape.ft, shape.dft)
    LU = lufact(A)

	ndft = sqrt.(sum(abs2,shape.dft,2))

    r = sqrt.((shape.ft[:,1] - ls_pos[1,1]).^2 + (shape.ft[:,2] - ls_pos[1,2]).^2)

    uls = -ls_amp[1]*0.25im*besselh.(0,kout*r)
    duls = ls_amp[1]*0.25im*(kout*besselh.(1,kout*r)).*((shape.ft[:,1]-ls_pos[1,1]).*shape.dft[:,2]-(shape.ft[:,2]-ls_pos[1,2]).*shape.dft[:,1])./r./ndft
    for i = 2:length(ls_amp)
        r = sqrt.((shape.ft[:,1] - ls_pos[i,1]).^2 + (shape.ft[:,2] - ls_pos[i,2]).^2)
        uls -= ls_amp[i]*0.25im*besselh.(0,kout*r)
        duls -= -ls_amp[i]*0.25im*kout*besselh.(1,kout*r).*((shape.ft[:,1]-ls_pos[i,1]).*shape.dft[:,2]-(shape.ft[:,2]-ls_pos[i,2]).*shape.dft[:,1])./r./ndft
    end

    #outer plane wave
    uinc = exp.(1.0im*kin*(cos(θ_i)*shape.ft[:,1] + sin(θ_i)*shape.ft[:,2]))
    duinc = (1.0im*kin)*(uinc.*(cos(θ_i)*shape.dft[:,2] - sin(θ_i)*shape.dft[:,1])./ndft)

    rhs = -[uinc+uls;duinc+duls]
    sigma_mu = LU\rhs
    return sigma_mu
end
