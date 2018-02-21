function plotNearField_pgf(filename, k0, kin, P, sp::ScatteringProblem, θ_i = 0;
                            opt::FMMoptions = FMMoptions(), use_multipole = true,
                            x_points = 201, y_points = 201, border = find_border(sp),
                            downsample = 1, include_preamble = false, normalize = 1.0)
    #important difference: here function is measured on the grid, unlike pcolormesh.
    shapes = sp.shapes;	ids = sp.ids; centers = sp.centers; φs = sp.φs
    x_min, x_max, y_min, y_max = border

    x = linspace(x_min, x_max, x_points)
    y = linspace(y_min, y_max, y_points)

    xgrid = repmat(x',y_points,1)
    ygrid = repmat(y,1,x_points)
    points = [xgrid[:] ygrid[:]]

    Ez = calc_near_field(k0, kin, P, sp, points, θ_i,
                            use_multipole=use_multipole, opt = opt)

    aspect = (y_max - y_min)/(x_max - x_min)

    dt = DataFrames.DataFrame(x = points[:,1]/normalize, y = points[:,2]/normalize, z = abs.(Ez))
    CSV.write(filename * ".dat", dt, delim = '\t')
    pgf.@pgf begin
        p = pgf.Plot3(pgf.Table(basename(filename) * ".dat"),
            {   surf,
                no_markers,
                shader = "interp",
                "mesh/rows" = y_points}, incremental = false)
        ax = pgf.Axis(p,
         {  xmin = x_min/normalize,
            xmax = x_max/normalize,
            ymin = y_min/normalize,
            ymax = y_max/normalize,
            xlabel = "\$x\$",
            ylabel = "\$y\$",
            scale_only_axis,
            width = "\\linewidth",
            height = "$aspect*\\linewidth",
            view = "{0}{90}",
            "mesh/ordering" = "y varies",
            colorbar,
            point_meta_max = maximum(abs.(Ez))})
    end
    drawShapes_pgf(shapes, centers, ids, φs, ax, 1.1*maximum(abs.(Ez)), downsample, normalize = normalize)
    pgf.save(filename, ax ,include_preamble = include_preamble)
end

function drawShapes_pgf(shapes, centers, ids, φs, ax, floating, downsample; normalize = 1.0)
    for ic = 1:size(centers,1)
        if typeof(shapes[ids[ic]]) == ShapeParams
            if φs[ic] == 0.0
                ft_rot = shapes[ids[ic]].ft[1:downsample:end,:] .+ centers[ic,:]'
            else
                Rot = [cos(φs[ic]) sin(φs[ic]);-sin(φs[ic]) cos(φs[ic])]
                ft_rot = shapes[ids[ic]].ft[1:downsample:end,:]*Rot .+ centers[ic,:]'
            end
            push!(ax, pgf.Plot3(pgf.Coordinates([ft_rot[:,1];ft_rot[1,1]]/normalize,
                [ft_rot[:,2];ft_rot[1,2]]/normalize, floating*ones(Float64,size(ft_rot,1)+1)),
                "black", "no markers", "thick", incremental = false))
        else
            x = centers[ic,1]/normalize
            y = centers[ic,2]/normalize
            R =  shapes[ids[ic]].R/normalize
            push!(ax, "\\addplot3 [black, thick, domain=0:2*pi,samples=100]
                        ({$x+$R*cos(deg(x))},{$y+$R*sin(deg(x))},$floating);")
        end
    end
end
