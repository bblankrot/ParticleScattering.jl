"""
    plot_near_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                        opt::FMMoptions = FMMoptions(), use_multipole = true,
                        x_points = 201, y_points = 201, border = find_border(sp),
                        normalize = 1.0)

Plots the total electric field as a result of a plane wave with incident
angle `θ_i` scattering from the ScatteringProblem `sp`, using matplotlib's
`pcolormesh`. Can accept number of sampling points in each direction plus
bounding box or calculate automatically.

Uses the FMM options given by `opt` (default behavious is disabled FMM);
`use_multipole` dictates whether electric field is calculated using the
multipole/cylindrical harmonics (true) or falls back on potential densities
(false). Either way, the multiple-scattering system is solved in the cylindrical
harmonics space. Normalizes all distances and sizes in plot (but not output) by
`normalize`.

Returns the calculated field in two formats:
1. `(points, Ez)` where `Ez[i]` is the total electric field at `points[i,:]`, and
2. `(xgrid,ygrid,zgrid)`, the format suitable for `pcolormesh`, where `zgrid[i,j]`
contains the field at `(mean(xgrid[i, j:j+1]), mean(ygrid[i:i+1, j]))`.
"""
function plot_near_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                        opt::FMMoptions = FMMoptions(), use_multipole = true,
                        x_points = 201, y_points = 201, interpolate = false,
                        border = find_border(sp), normalize = 1.0)

    x_min, x_max, y_min, y_max = border

    # with pcolormesh, the field is calculated at the center of each rectangle.
    # thus we need two grids - the rectangle grid and the sampling grid.
    x = linspace(x_min, x_max, x_points + 1)
    y = linspace(y_min, y_max, y_points + 1)

    xgrid = repmat(x', y_points + 1, 1)
    ygrid = repmat(y, 1, x_points + 1)
    dx = (x_max - x_min)/2/x_points
    dy = (y_max - y_min)/2/y_points
    points = cat(2, vec(xgrid[1:y_points, 1:x_points]) + dx,
                vec(ygrid[1:y_points, 1:x_points]) + dy)

    if interpolate == true
        #TODO: make this work, some time
        flags = Array{Bool}(size(points,1))
        flags[:] = true
        for ix = 1:size(points,1)
            for ic = 1:length(sp.ids)
                typ = sqrt(sum(abs2,sp.shapes[sp.ids[ic]].ft[1,:] - sp.shapes[sp.ids[ic]].ft[2,:]))
                dist = minimum(sqrt(sum(abs2,(points[ix,:] - sp.centers[ic,:])' .- sp.shapes[sp.ids[ic]].ft ,2)))
                if dist < 3*typ
                    flags[ix] = false
                    break
                end
            end
        end
        points = points[flags,:]
        Ez = calc_near_field(k0, kin, P, sp, points, θ_i, use_multipole=use_multipole, opt = opt)
        levs = linspace(0.0,maximum(abs.(Ez)),40)
        x_min += dx; y_min += dy; x_max += dx; y_max += dy;
        figure()
        tricontourf(vec(points[:,1]),vec(points[:,2]),vec(abs.(Ez)),levels=levs)
    else
        Ez = calc_near_field(k0, kin, P, sp, points, θ_i,
                use_multipole=use_multipole, opt = opt)
        zgrid = reshape(Ez, y_points, x_points)
        figure()
        pcolormesh(xgrid/normalize, ygrid/normalize, abs.(zgrid))
    end

    ax = gca()
    draw_shapes(sp.shapes, sp.centers, sp.ids, sp.φs, ax; normalize = normalize)
    xlim([x_min/normalize;x_max/normalize])
    ylim([y_min/normalize;y_max/normalize])
    tight_layout()
    ax[:set_aspect]("equal", adjustable = "box")
    return (points,Ez),(xgrid,ygrid,zgrid)
end

"""
    plot_far_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                        opt::FMMoptions = FMMoptions(), use_multipole = true,
                        plot_points = 200)

Plots the echo width (radar cross section in two dimensions) for a given
scattering problem. `opt`, `use_multipole` are as in `plot_near_field`. Also
returns the echo width.
"""
function plot_far_field(k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                    opt::FMMoptions = FMMoptions(), use_multipole = true,
                    plot_points = 200)

    Rmax = maximum(s.R for s in sp.shapes)

    x_max,y_max = maximum(sp.centers,1) + 2*Rmax
    x_min,y_min = minimum(sp.centers,1) - 2*Rmax
    Raggregate = 0.5*max(x_max - x_min, y_max - y_min)
    x_center = 0.5*(x_max + x_min)
    y_center = 0.5*(y_max + y_min)
    Rfar = Raggregate*1e6
    theta_far = linspace(0, 2π, plot_points)
    x_far = x_center + Rfar*cos.(theta_far)
    y_far = y_center + Rfar*sin.(theta_far)
    points = [x_far y_far]

    Ez = calculateFarField(k0, kin, P, points, sp, θ_i,
            use_multipole = use_multipole, opt = opt)
    Ez[:] = (k0*Rfar)*abs2.(Ez)
    #plot echo width
    figure()
    plot(theta_far/π, Ez)
    xlabel("\$ \\theta/\\pi \$")
    ylabel("\$ \\sigma/\\lambda_0 \$")
    title("Echo Width")
    tight_layout()
    xlim([0;2])
    return Ez
end

"""
    draw_shapes(shapes, centers, ids, φs, ax = gca())

Draws all of the shapes in a given scattering problem. Parametrized shapes are
drawn as polygons while circles are drawn using matplotlib's `patch.Circle`.
"""
function draw_shapes(shapes, centers, ids, φs, ax = gca(); normalize = 1.0)
    #draw shapes
    for ic = 1:size(centers,1)
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = shapes[ids[ic]].ft .+ centers[ic,:]'
            else
                Rot = cartesianrotation(φs[ic])
                ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
            end
            ax[:plot]([ft_rot[:,1];ft_rot[1,1]]/normalize,
                    [ft_rot[:,2];ft_rot[1,2]]/normalize, "k", linewidth = 2)
        else
            ax[:add_patch](patch.Circle((centers[ic,1]/normalize,
                            centers[ic,2]/normalize),
                            radius = shapes[ids[ic]].R/normalize,
                            edgecolor="k", facecolor="none", linewidth = 2))
        end
    end
end

"""
    calc_near_field(k0, kin, P, sp::ScatteringProblem, points, θ_i;
                            opt::FMMoptions = FMMoptions(), use_multipole = true,
                            verbose = true)

Calculates the total electric field as a result of a plane wave with incident
angle `θ_i` scattering from the ScatteringProblem `sp`, at `points`.
Uses the FMM options given by `opt` (default behavious is disabled FMM);
`use_multipole` dictates whether electric field is calculated using the
multipole/cylindrical harmonics (true) or falls back on potential densities
(false). Either way, the multiple-scattering system is solved in the cylindrical
harmonics space, and the field by a particular scatterer inside its own scattering
discs is calculated by potential densities, as the cylindrical harmonics
approximation is not valid there.
"""
function calc_near_field(k0, kin, P, sp::ScatteringProblem, points, θ_i;
                            opt::FMMoptions = FMMoptions(), use_multipole = true,
                            verbose = true)

    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    u = zeros(Complex{Float64},size(points,1))
    if opt.FMM
        result,sigma_mu =  solve_particle_scattering_FMM(k0, kin, P, sp, θ_i, opt,
                            verbose = verbose)
        if result[2].isconverged == false
            warn("FMM process did not converge")
            return
        end
        beta = result[1]
    else
        beta, sigma_mu = solve_particle_scattering(k0, kin, P, sp, θ_i)
    end

    tic()
    #first, let's mark which points are in which shapes in tags:
    #0 denotes outside everything, +-i means inside shape i or between it and its "multipole disk"
    tags = tagpoints(sp, points)
    dt_tag = toq()

    tic()
    rng_in = zeros(Bool,size(points,1))
    rng_out = zeros(Bool,size(points,1))
    Rot = Array{Float64}(2,2)
    for ic = 1:size(sp)
        rng_in[:] = (tags .== ic)
        rng_out[:] = (tags .== -ic)
        (any(rng_in) || any(rng_out)) || continue
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = shapes[ids[ic]].ft .+ centers[ic,:]'
                dft_rot = shapes[ids[ic]].dft
            else
                Rot[:] = cartesianrotation(φs[ic])
                ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                dft_rot = shapes[ids[ic]].dft*Rot.'
            end
        end
        #field inside shape
        if any(rng_in)
            if typeof(shapes[ids[ic]]) == ShapeParams
                u[rng_in] += scatteredField(sigma_mu[ic], kin, shapes[ids[ic]].t,
                                ft_rot, dft_rot, points[rng_in,:])
            else
                u[rng_in] += innerFieldCircle(kin, sigma_mu[ic], centers[ic,:],
                                points[rng_in,:])
            end
        end
        #field between shape and multipole disk (impossible for circle)
        if any(rng_out)
            u[rng_out] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t,
                            ft_rot, dft_rot, points[rng_out,:])
            u[rng_out] += exp.(1.0im*k0*(cos(θ_i)*points[rng_out,1] +
                                        sin(θ_i)*points[rng_out,2])) #incident
            for ic2 = 1:size(sp)
                ic == ic2 && continue
                if use_multipole
                    scattered_field_multipole!(u, k0, beta, P, centers, ic2, points,
                        find(rng_out))
                else
                    if φs[ic2] == 0.0
                        ft_rot2 = shapes[ids[ic2]].ft .+ centers[ic2,:]'
                        dft_rot2 = shapes[ids[ic2]].dft
                    else
                        Rot[:] = cartesianrotation(φs[ic2])
                        ft_rot2 = shapes[ids[ic2]].ft*Rot.' .+ centers[ic2,:]'
                        dft_rot2 = shapes[ids[ic2]].dft*Rot.'
                    end
                    u[rng_out] += scatteredField(sigma_mu[ic2], k0, shapes[ids[ic2]].t,
                                    ft_rot2, dft_rot2, points[rng_out,:])
                end
            end
        end
    end
    dt_in = toq()

    tic()
    #now compute field outside all shapes
    rng = (tags .== 0)
    #incident field
    u[rng] = exp.(1.0im*k0*(cos(θ_i)*points[rng,1] + sin(θ_i)*points[rng,2]))
    if use_multipole
        scattered_field_multipole!(u, k0, beta, P, centers, 1:size(sp), points, find(rng))
    else
        for ic = 1:size(centers,1)
            if typeof(shapes[ids[ic]]) == ShapeParams
                if φs[ic] == 0.0
                    ft_rot = shapes[ids[ic]].ft .+ centers[ic,:]'
                    u[rng] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t,
                                ft_rot, shapes[ids[ic]].dft, points[rng,:])
                else
                    Rot[:] = cartesianrotation(φs[ic])
                    ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                    dft_rot = shapes[ids[ic]].dft*Rot.'
                    u[rng] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t,
                                ft_rot, dft_rot, points[rng,:])
                end
            else
                scattered_field_multipole!(u, k0, beta, P, centers, ic, points,
                    find(rng))
            end
        end
    end
    dt_out = toq()
    if verbose
        println("Time spent calculating field:")
        println("Location tagging: $dt_tag")
        println("In/around scatterers: $dt_in")
        println("Outside scatterers: $dt_out")
    end
    return u
end

function calculateFarField(k0, kin, P, points, sp::ScatteringProblem, θ_i;
                        opt::FMMoptions = FMMoptions(), use_multipole = true)
    #calc only scattered field + assumes all points are outside shapes
    shapes = sp.shapes; centers = sp.centers; ids = sp.ids, φs = sp.φs
    if opt.FMM
        result,sigma_mu =  solve_particle_scattering_FMM(k0, kin, P, sp, θ_i, opt)
        if result[2].isconverged == false
            warn("FMM process did not converge")
            return
        end
        beta = result[1]
    else
        beta, sigma_mu = solve_particle_scattering(k0, kin, P, sp, θ_i)
    end
    Ez = zeros(Complex{Float64}, size(points,1))
    if use_multipole
        scattered_field_multipole!(Ez, k0, beta, P, centers, 1:size(sp),
            points, 1:size(points,1))
    else
        for ic = 1:size(sp)
            if typeof(shapes[ids[ic]]) == ShapeParams
                if φs[ic] == 0.0
                    ft_rot = shapes[ids[ic]].ft .+ centers[ic,:]'
                    Ez[:] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t,
                                ft_rot, shapes[ids[ic]].dft, points)
                else
                    Rot = cartesianrotation(φs[ic])
                    ft_rot = shapes[ids[ic]].ft*Rot.' .+ centers[ic,:]'
                    dft_rot = shapes[ids[ic]].dft*Rot.'
                    Ez[:] += scatteredField(sigma_mu[ic], k0, shapes[ids[ic]].t,
                                ft_rot, dft_rot, points)
                end
            else
                warning("should be only subset of beta!")
                scattered_field_multipole!(Ez, k0, beta, P, centers, ic, points,
                    1:size(points,1))
            end
        end
    end
    return Ez
end

function tagpoints(sp, points)
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs

    tags = zeros(Integer, size(points,1))
    X = Array{Float64}(2)
    for ix = 1:size(points,1)
        for ic = 1:size(sp)
            X .= points[ix,:] - centers[ic,:]
            if sum(abs2,X) <= shapes[ids[ic]].R^2
                if typeof(shapes[ids[ic]]) == ShapeParams
                    if φs[ic] != 0.0 #rotate point backwards instead of shape forwards
                        Rot = [cos(-φs[ic]) -sin(-φs[ic]);sin(-φs[ic]) cos(-φs[ic])]
                        X = Rot*X
                    end
                    tags[ix] = pInPolygon(X, shapes[ids[ic]].ft) ? ic : -ic
        			break #can't be in two shapes
                else #CircleParams
                    tags[ix] = ic
                    break #can't be in two shapes
                end
        	end
        end
    end
    tags
end
