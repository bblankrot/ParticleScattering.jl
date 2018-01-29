function plotNearField_pgf(filename, k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                            opt::FMMoptions = FMMoptions(), use_multipole = true,
                            x_points = 201, y_points = 201, border = find_border(sp),
                            downsample = 1, include_preamble = false)
    #important difference: here function is measured on the grid, unlike pcolormesh.
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    axwidth = 5
    
    x_min, x_max, y_min, y_max = border

    x = linspace(x_min, x_max, x_points)
    y = linspace(y_min, y_max, y_points)

    xgrid = repmat(x',y_points,1)
    ygrid = repmat(y,1,x_points)
    points = [xgrid[:] ygrid[:]]

    Ez = calc_near_field(k0, kin, P, sp, points, θ_i,
                            use_multipole=use_multipole, opt = opt)

    normalized = true
    l0 = (normalized ? 2π/k0 : 1.0)
    aspect = (y_max - y_min)/(x_max - x_min)

    dt = DataFrames.DataFrame(x = points[:,1], y = points[:,2], z = abs.(Ez))
    CSV.write(filename * ".dat", dt, delim = '\t')
    pgf.@pgf begin
        p = pgf.Plot3(pgf.Table(basename(filename) * ".dat"),
            {   surf,
                no_markers,
                shader = "interp",
                "mesh/rows" = y_points}, incremental = false)
        if normalized
            ax = pgf.Axis(p,
             {  xmin = x_min/l0,
                xmax = x_max/l0,
                ymin = y_min/l0,
                ymax = y_max/l0,
                xlabel = "\$x/\\lambda_0\$",
                ylabel = "\$y/\\lambda_0\$",
                scale_only_axis,
                width = "\\figurewidth",
                height = "$aspect*\\figurewidth",
                view = "{0}{90}",
                "mesh/ordering" = "y varies",
                colorbar,
                point_meta_max = maximum(abs.(Ez))})
        else
            ax = pgf.Axis(p,
             {  xmin = x_min,
                xmax = x_max,
                ymin = y_min,
                ymax = y_max,
                xlabel = "\$x\$",
                ylabel = "\$y\$",
                scale_only_axis,
                width = "\\figurewidth",
                height = "$aspect*\\figurewidth",
                view = "{0}{90}",
                "mesh/ordering" = "y varies",
                colorbar,
                point_meta_max = maximum(abs.(Ez))})
        end
    end
    drawShapes_pgf(l0, shapes, centers, ids, φs, ax, 1.1*maximum(abs.(Ez)), downsample)
    pgf.save(filename, ax ,include_preamble = include_preamble)
end

function drawShapes_pgf(l0, shapes, centers, ids, φs, ax, floating, downsample)
    for ic = 1:size(centers,1)
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = shapes[ids[ic]].ft[1:downsample:end,:] .+ centers[ic,:]'
            else
                Rot = [cos(φs[ic]) sin(φs[ic]);-sin(φs[ic]) cos(φs[ic])]
                ft_rot = shapes[ids[ic]].ft[1:downsample:end,:]*Rot .+ centers[ic,:]'
            end
            push!(ax, pgf.Plot3(pgf.Coordinates([ft_rot[:,1];ft_rot[1,1]]/l0,
                [ft_rot[:,2];ft_rot[1,2]]/l0, floating*ones(Float64,size(ft_rot,1)+1)),
                "black", "no markers", "thick", incremental = false))
        else
            x = centers[ic,1]/l0
            y = centers[ic,2]/l0
            R =  shapes[ids[ic]].R/l0
            push!(ax, "\\addplot3 [black, thick, domain=0:2*pi,samples=100]
                        ({$x+$R*cos(deg(x))},{$y+$R*sin(deg(x))},$floating);")
        end
    end
end
