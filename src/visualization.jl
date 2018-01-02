function plotNearField(k0, kin, P, sp::ScatteringProblem, θ_i = 0; opt::FMMoptions = FMMoptions(), use_multipole = true, x_points = 201, y_points = 201, interpolate = false, border = [NaN64;NaN64;NaN64;NaN64])
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    if any(isnan.(border))
        Rmax = 0.0
        for s in shapes
            Rs = (typeof(s) == ShapeParams ? maximum(hypot.(s.ft[:,1],s.ft[:,2])) : s.R)
            Rmax = max(Rs,Rmax)
        end
        (x_max,y_max) = maximum(centers,1) + 2*Rmax
        (x_min,y_min) = minimum(centers,1) - 2*Rmax
    else
        x_min = border[1]; x_max = border[2]
        y_min = border[3]; y_max = border[4]
    end

    x = linspace(x_min, x_max, x_points + 1)
    y = linspace(y_min,y_max,y_points + 1)

    xgrid = repmat(x',y_points + 1,1)
    ygrid = repmat(y,1,x_points + 1)
    dx = (x_max - x_min)/2/x_points
    dy = (y_max - y_min)/2/y_points
    #so the color of each square is the the value in its center
    xgrid2 = repmat(x[1:x_points]',y_points,1)
    ygrid2 = repmat(y[1:y_points],1,x_points)
    points = [xgrid2[:]+dx ygrid2[:]+dy]

    if interpolate == true
        #TODO: make this work, some time
        flags = Array{Bool}(size(points,1))
        flags[:] = true
        for ix = 1:size(points,1)
            for ic = 1:length(ids)
                typ = sqrt(sum(abs2,shapes[ids[ic]].ft[1,:] - shapes[ids[ic]].ft[2,:]))
                dist = minimum(sqrt(sum(abs2,(points[ix,:] - centers[ic,:])' .- shapes[ids[ic]].ft ,2)))
                if dist < 3*typ
                    flags[ix] = false
                    break
                end
            end
        end
        points = points[flags,:]
        Ez = calculateNearField(k0, kin, P, sp, points, θ_i, use_multipole=use_multipole, opt = opt)
        levs = linspace(0.0,maximum(abs.(Ez)),40)
        x_min += dx; y_min += dy; x_max += dx; y_max += dy;
        figure()
        tricontourf(vec(points[:,1]),vec(points[:,2]),vec(abs.(Ez)),levels=levs)
    else
        Ez = calculateNearField(k0, kin, P, sp, points, θ_i, use_multipole=use_multipole, opt = opt)
        zgrid = reshape(Ez, y_points, x_points)
        #plot pcolormesh of absolute Ez
        figure()
        pcolormesh(xgrid,ygrid,abs.(zgrid))
    end

    ax = gca()
    drawShapes(shapes, centers, ids, φs, ax)
    # xlabel(L"x")
    # ylabel(L"y")
    # title(L"|E_z|")
    xlim([x_min;x_max])
    ylim([y_min;y_max])
    # colorbar()
    tight_layout()
    ax[:set_aspect]("equal", adjustable = "box")
    return points,Ez,(xgrid,ygrid,zgrid)
end

function plotFarField(k0, kin, P, sp::ScatteringProblem, θ_i = 0; opt::FMMoptions = FMMoptions(), use_multipole = true, plot_points = 200)
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    Rmax = 0.0
    for s in shapes
        Rs = (typeof(s) == ShapeParams ? maximum(hypot.(s.ft[:,1],s.ft[:,2])) : s.R)
        Rmax = max(Rs,Rmax)
    end
    x_min = minimum(centers[:,1]) - Rmax;
    x_max = maximum(centers[:,1]) + Rmax;
    y_min = minimum(centers[:,2]) - Rmax;
    y_max = maximum(centers[:,2]) + Rmax;
    Raggregate = 0.5*max(x_max - x_min, y_max - y_min)
    x_center = 0.5*(x_max + x_min)
    y_center = 0.5*(y_max + y_min)
    Rfar = Raggregate*1e6
    theta_far = linspace(0,2*pi,plot_points)
    x_far = x_center + Rfar*cos.(theta_far)
    y_far = y_center + Rfar*sin.(theta_far)
    points = [x_far y_far]

    Ez = calculateFarField(k0, kin, P, points, centers, φs, ids, shapes, θ_i, use_multipole=use_multipole, opt = opt)

    #plot rcs
    figure()
    plot(theta_far/pi,2*pi*Rfar*abs.(Ez).^2/(2*pi/k0))
    xlabel(L"\theta/\pi")
    ylabel(L"\sigma/\lambda_0")
    title("Echo Width")
    tight_layout()
    xlim([0;2])
    return Ez
end

function drawShapes(shapes, centers, ids, φs, ax = gca())
    #draw shapes
    for ic = 1:size(centers,1)
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = shapes[ids[ic]].ft .+ centers[ic,:]'
            else
                Rot = [cos(φs[ic]) sin(φs[ic]);-sin(φs[ic]) cos(φs[ic])]
                ft_rot = shapes[ids[ic]].ft*Rot .+ centers[ic,:]'
            end
            ax[:plot]([ft_rot[:,1];ft_rot[1,1]],[ft_rot[:,2];ft_rot[1,2]],
                        "k", linewidth = 2)
        else
            ax[:add_patch](patch.Circle((centers[ic,1],centers[ic,2]),
                            radius = shapes[ids[ic]].R, edgecolor="k",
                            facecolor="none", linewidth = 2))
        end
    end
end

function calculateNearField(k0, kin, P, sp::ScatteringProblem, points, θ_i; opt::FMMoptions = FMMoptions(), use_multipole = true)
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    z = zeros(Complex{Float64},size(points,1))
    if opt.FMM
        result,sigma_mu =  solveParticleScattering_FMM(k0, kin, P, sp, θ_i, opt)
        if result[2].isconverged == false
            warn("FMM process did not converge")
            return
        end
        beta = result[1]
    else
        (beta, sigma_mu) = solveParticleScattering(k0, kin, P, sp, θ_i)
    end
    tic()
    #first, let's mark which points are in which shapes in tags:
    #0 denotes outside everything, +-i means inside shape i or between it and its "multipole disk"
    rng_in = zeros(Bool,size(points,1))
    rng_out = zeros(Bool,size(points,1))
    tags = zeros(Int64, size(points,1))
    for ix = 1:size(points,1)
    	for ic = 1:size(centers,1)
            X2 = points[ix,:] - centers[ic,:]
            if sum(abs2,X2) <= shapes[ids[ic]].R^2
                if typeof(shapes[ids[ic]]) == ShapeParams
                    if φs[ic] != 0.0 #rotate backwards
                        Rot = [cos(-φs[ic]) -sin(-φs[ic]);sin(-φs[ic]) cos(-φs[ic])]
                        X2 = Rot*X2
                    end
                    tags[ix] = pInPolygon(X2, shapes[ids[ic]].ft) ? ic : -ic
        			break #can't be in two shapes
                else #CircleParams
                    tags[ix] = ic
                    break #can't be in two shapes
                end
    		end
    	end
    end
    dt_tag = toq()
    tic()
    for ic = 1:size(centers,1)
        rng_in[:] = (tags .== ic)
        rng_out[:] = (tags .== -ic)
        (any(rng_in) || any(rng_out)) || continue
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = [shapes[ids[ic]].ft[ii,jj] + centers[ic,jj] for ii=1:size(shapes[ids[ic]].ft,1), jj=1:2]
                dft_rot = shapes[ids[ic]].dft
            else
                Rot = [cos(φs[ic]) -sin(φs[ic]);sin(φs[ic]) cos(φs[ic])]
                ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                dft_rot = shapes[ids[ic]].dft*Rot.'
            end
        end
        #field inside shape
        if any(rng_in)
            if typeof(shapes[ids[ic]]) == ShapeParams
                z[rng_in] += scatteredField(sigma_mu[ic], kin, shapes[ids[ic]].t, ft_rot, dft_rot, points[rng_in,:])
            else
                z[rng_in] += innerFieldCircle(kin, sigma_mu[ic], centers[ic,:], points[rng_in,:])
            end
        end
        #field between shape and multipole disk (impossible for circle)
        if any(rng_out)
            z[rng_out] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t, ft_rot, dft_rot, points[rng_out,:])
            z[rng_out] += exp.(1.0im*k0*(cos(θ_i)*points[rng_out,1] + sin(θ_i)*points[rng_out,2])) #incident
            for ic2 = 1:size(centers,1)
                ic == ic2 && continue
                if use_multipole
                    scatteredFieldMultipole(k0, beta, P, centers, ic2, points, z, find(rng_out))
                else
                    if φs[ic2] == 0.0
                        ft_rot2 = [shapes[ids[ic2]].ft[ii,jj] + centers[ic2,jj] for ii=1:size(shapes[ids[ic2]].ft,1), jj=1:2]
                        dft_rot2 = shapes[ids[ic2]].dft
                    else
                        Rot = [cos(φs[ic2]) -sin(φs[ic2]);sin(φs[ic2]) cos(φs[ic2])]
                        ft_rot2 = shapes[ids[ic2]].ft*Rot.' .+ centers[ic2,:]'
                        dft_rot2 = shapes[ids[ic2]].dft*Rot.'
                    end
                    z[rng_out] += scatteredField(sigma_mu[ic2], k0, shapes[ids[ic2]].t, ft_rot2, dft_rot2, points[rng_out,:])
                end
            end
        end
    end
    dt_in = toq()
    tic()
    #now compute field outside all shapes
    rng = (tags .== 0)
    einc = exp.(1.0im*k0*(cos(θ_i)*points[rng,1] + sin(θ_i)*points[rng,2]))
    z[rng] = einc
    if use_multipole
        scatteredFieldMultipole(k0, beta, P, centers, 1:size(centers,1), points, z, find(rng))
    else
        for ic = 1:size(centers,1)
            if typeof(shapes[ids[ic]]) == ShapeParams
                if φs[ic] == 0.0
                    ft_rot = [shapes[ids[ic]].ft[ii,jj] + centers[ic,jj] for ii=1:size(shapes[ids[ic]].ft,1), jj=1:2]
                    z[rng] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t, ft_rot, shapes[ids[ic]].dft, points[rng,:])
                else
                    Rot = [cos(φs[ic]) -sin(φs[ic]);sin(φs[ic]) cos(φs[ic])]
                    ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                    dft_rot = shapes[ids[ic]].dft*Rot.'
                    z[rng] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t, ft_rot, dft_rot, points[rng,:])
                end
            else
                warning("should be only subset of beta!")
                scatteredFieldMultipole(k0, beta, P, centers, ic, points, z, find(rng))
            end
        end
    end
    dt_out = toq()
    println("Time spent calculating field:")
    println("Location tagging: $dt_tag")
    println("In/around scatterers: $dt_in")
    println("Outside scatterers: $dt_out")
    return z
end

function calculateFarField(k0, kin, P, points, centers, φs, ids, shapes, θ_i; opt::FMMoptions = FMMoptions(), use_multipole = true)
    #calc only scattered field + assumes all points are outside shapes
    sp = ScatteringProblem(shapes,ids,centers,φs)
    if opt.FMM
        result,sigma_mu =  solveParticleScattering_FMM(k0, kin, P, sp, θ_i, opt)
        if result[2].isconverged == false
            warn("FMM process did not converge")
            return
        end
        beta = result[1]
    else
        (beta, sigma_mu) = solveParticleScattering(k0, kin, P,sp, θ_i)
    end
    Ez = zeros(Complex{Float64},size(points,1))
    if use_multipole
        scatteredFieldMultipole(k0, beta, P, centers, 1:size(centers,1), points, Ez, 1:size(points,1))
    else
        for ic = 1:size(centers,1)
            if typeof(shapes[ids[ic]]) == ShapeParams
                if φs[ic] == 0.0
                    ft_rot = [shapes[ids[ic]].ft[ii,jj] + centers[ic,jj] for ii=1:size(shapes[ids[ic]].ft,1), jj=1:2]
                    Ez[:] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t, ft_rot, shapes[ids[ic]].dft, points)
                else
                    Rot = [cos(φs[ic]) -sin(φs[ic]);sin(φs[ic]) cos(φs[ic])]
                    ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                    dft_rot = shapes[ids[ic]].dft*Rot.'
                    Ez[:] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t, ft_rot, dft_rot, points)
                end
            else
                warning("should be only subset of beta!")
                scatteredFieldMultipole(k0, beta, P, centers, ic, points, Ez, 1:size(points,1))
            end
        end
    end
    return Ez
end
