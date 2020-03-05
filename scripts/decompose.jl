using HDF5
using MapImages
using MAT

function solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
    # @show size(asc_los_map), size(desc_los_map), size(asc_img), size(desc_img)

    east = similar(asc_img)
    up = similar(asc_img)

    for jj = 1:size(asc_los_map, 2)
        for ii = 1:size(asc_los_map, 1)
            asc_eu = asc_los_map[ii, jj, [1, 3]]
            desc_eu = desc_los_map[ii, jj, [1, 3]]

            A = hcat(asc_eu, desc_eu)'
            b = [asc_img[ii, jj] ; desc_img[ii, jj]]

            b = asc_eu[2] < 0 ? -b : b

            x = A \ b
            east[ii, jj] = x[1]
            up[ii, jj] = x[2]
        end
    end
    return east, up
end

LOS_MAP = "los_map.h5"  

function solve_east_up(asc_path::AbstractString, desc_path::AbstractString,
                       asc_fname::AbstractString, desc_fname::AbstractString=asc_fname, 
                       dset="velos/1")
    asc_img = permutedims(h5read(joinpath(asc_path, asc_fname), dset))
    desc_img = permutedims(h5read(joinpath(desc_path, desc_fname), dset))

    asc_img, desc_img = MapImages._mask_asc_desc(asc_img, desc_img)

    asc_los_map = permutedims(h5read(joinpath(asc_path, LOS_MAP), "stack"), (2, 1, 3))
    desc_los_map =permutedims(h5read(joinpath(desc_path, LOS_MAP), "stack"), (2, 1, 3))
    return solve_east_up(asc_img, desc_img, asc_los_map, desc_los_map)
end

function plot_east_up(east, up; cmap="seismic_wide_y", vm=20, east_scale=1.0, title="", show_cbar=false)
    fig, axes = plt.subplots(1, 2, sharex=true, sharey=true)
    axim1 = axes[1].imshow(up, cmap=cmap, vmin=-vm, vmax=vm)
    # fig.colorbar(axim1, ax=axes[1])

    vmeast = east_scale * vm
    axim2 = axes[2].imshow(east, cmap=cmap, vmin=-vmeast, vmax=vmeast)
    show_cbar && fig.colorbar(axim2, ax=axes[2])

    axes[1].set_title("up")
    axes[2].set_title("east")
    fig.suptitle(title)
    return fig, axes
end

function demo_east_up(fn="velocities_prune_l1.h5"; dset="velos/1", full=true, vm=20, east_scale=1.0, shifta=0.0, shiftd=0.0, 
                      asc_path="/data1/scott/pecos/path78-bbox2/igrams_looked_18/", desc_path="/data4/scott/path85/stitched/igrams_looked_18/",
                      cmap="seismic_wide_y", show=true)
    if full
        # asc_path, desc_path = ("/data1/scott/pecos/path78-bbox2/igrams_looked/", "/data4/scott/path85/stitched/igrams_looked/")
        asc_fname, desc_fname = map(x -> joinpath(x, fn), (asc_path, desc_path))
        asc_img, desc_img = MapImages.find_overlaps(asc_fname, desc_fname, dset)

        asc_los_map = MapImage(joinpath(asc_path, LOS_MAP), dset_name="stack")
        desc_los_map = MapImage(joinpath(desc_path, LOS_MAP), dset_name="stack")
        #
        asc_idxs, desc_idxs = MapImages.find_overlap_idxs(asc_los_map, desc_los_map)
        # @show size(asc_los_map[asc_idxs..., :]), size(desc_los_map[desc_idxs..., :])

        east, up = solve_east_up(asc_img .+ shifta, desc_img .+ shiftd, asc_los_map[asc_idxs..., :], desc_los_map[desc_idxs..., :])
    else
        asc_path, desc_path = ("/data3/scott/pecos/zoom_pecos_full_78/igrams_looked/", "/data3/scott/pecos/zoom_pecos_full_85/igrams_looked/", )
        east, up = solve_east_up(asc_path, desc_path, fn, fn, dset)
    end
    if show
        fig, axes = plot_east_up(east, up; cmap=cmap, vm=vm, title="$fn: $dset", east_scale=east_scale)
    end
    # return east, up, fig, axes
    return east, up
end


function eastup_insar(stations, east, up)
    easts, ups = [], []
    for s in stations
        rc = gps.station_rowcol(s, Dict(up.demrsc))
        push!(easts, east[rc...])
        push!(ups, up[rc...])
    end
    return easts, ups
end

function slope_mm_yr(dts, data)
    # Convert to "days since start" for line fitting
    gps_poly = fit_line(dts, data)
    slope = length(gps_poly) == 2 ? Polynomials.coeffs(gps_poly)[2] : Polynomials.coeffs(gps_poly)[1]
    return 365 * 10 * slope
end

function eastup_gps(stations, start_date=Date(2014,11,1), end_date=Date(2019,1,1))

    easts, ups = [], []
    for s in stations
        dts, east, north, up = get_gps_enu(s)
        push!(easts, slope_mm_yr(dts, east))
        push!(ups, slope_mm_yr(dts, up))
    end
    return easts, ups
end

function eastup_diffs(stations, east, up, start_date=Date(2014,11,1), end_date=Date(2019,1,1))
    return eastup_insar(stations, east, up) .- eastup_gps(stations, start_date, end_date)
end

function save_east_up_decomp(fname, shifta, shiftd; vm=12, show=true, dset="velos/1", 
                             outname="zoom_pecos_vertical_east.mat",
                             latrange = (30.901666673336, 31.90000000028),
                             lonrange = (-104.0, -103.001666673056))
# lats = read(ff, "lats");
    # lons = read(ff, "lons");
    east, up = demo_east_up(fname; full=true, dset=dset,  shifta=shifta, shiftd=shiftd, vm=vm, show=show)
    g = Sario.load_geolist_from_h5(fname, dset)
    days = (g[end] - g[1]).value
    @show days

    if isnothing(latrange) || isnothing(lonrange)
        eastcut = east .* (days  / 3650)
        upcut = up .* (days  / 3650)
    else
        eastcut = east[latrange, lonrange] .* (days  / 3650)
        upcut = up[latrange, lonrange] .* (days  / 3650)
    end
    gx, gy = MapImages.grid(eastcut.demrsc)

    println("Lat range: $(extrema(gy))")
    println("Lon range: $(extrema(gx))")
    println("Saving lats, lons, ups, easts to $outname")

    dd = Dict("lats"=> gy, "lons" => gx, 
                           "up" => upcut.image, 
                           "east" => eastcut.image)
    isnothing(outname) && return dd
    matwrite(outname, dd);
    return dd
end

function save_all_east_up(outname=outname)
    dd_all = Dict{String, Any}()

    # p78 = "/data1/scott/pecos/path78-bbox2/igrams_looked/"
    # p85 = "/data4/scott/path85/stitched/igrams_looked/"

    varnames, values = [], []
    for year in 2016:2018
        fname = "velocities_$(year)_current.h5"
        dd = save_east_up_decomp(fname, 0, 0, dset="velos_shifted/1", show=false, outname=nothing, latrange=nothing)
        dd_all["lats"] = dd["lats"]
        dd_all["lons"] = dd["lons"]
        dd_all["east_$year"] = dd["east"]
        dd_all["up_$year"] = dd["up"]
        dd_all["east_$year"] = dd["east"]
    end
    matwrite(outname , dd_all);
    return dd_all
end

station_overlap = ["TXMH", "TXFS", "TXAD", "TXS3", "NMHB", "TXKM"]
eeups(fname, shifta, shiftd, dset="velos/1", full=true; show=true) = extrema.(eastup_insar(station_overlap, 
                                                                           demo_east_up(fname ;full=full, dset=dset, shifta=shifta, shiftd=shiftd, show=show)...))
eeups_diff(fname, shifta, shiftd, dset="velos/1", full=true; show=true) = eastup_diffs(station_overlap, demo_east_up(fname ;full=full, dset=dset, shifta=shifta, shiftd=shiftd, show=show)...)

# h5write("zoom_pecos_vertical_east.h5", "east_18", permutedims(read(ff3, "east")))
