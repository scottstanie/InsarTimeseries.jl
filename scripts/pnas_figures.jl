
# Formatting functions:
_remove_ticks(ax) = ax.tick_params(axis="both", which="both", bottom=false, top=false, labelbottom=false,
                                   left=false, right=false, labelleft=false)
_set_figsize(fig, h=h, w=w) = (fig.set_figheight(h), fig.set_figwidth(w))

function pnas_outlier_figure_old(station_name="TXMC", outlier_color="r", h=3, w=3.5, bins=40)

    # unw_vals = [h5read("gps_pixels.h5", name) for name in station_name_list78];
    unw_vals = h5read("gps_pixels_78.h5", station_name)
    geolist = Sario.load_geolist_from_h5("gps_pixels_78.h5");
    intlist = Sario.load_intlist_from_h5("gps_pixels_78.h5");

    end_date = Date(2018,1,1)
    geolist18 = geolist[geolist .> end_date]
    idxs = InsarTimeseries._good_idxs(geolist18, intlist)
    intlist17 = intlist[idxs]

    geolist17 = geolist[geolist .< end_date]
    unw_vals17 = unw_vals[idxs]
    
    # 1. plot grouped by date
    fig, ax = plot_big_days(geolist17, intlist17, unw_vals17, legend=false, color=outlier_color)
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers2.pdf", bbox_inches="tight", transparent=true, dpi=100)

    # 2. cm vs baseline igram outlier plot
    # function plot_grouped_by_day(geo, int, vals, nsigma=0, min_spread=2)
    geo = geolist17
    means = InsarTimeseries.mean_abs_val(geo, intlist17, p2c * unw_vals17)
    println("median val: $(median(means))")
    fig, ax = plt.subplots()
    ax.scatter(geo, means, label="means by date")
    nsigma = 4
    low, high = InsarTimeseries.two_way_cutoff(means, nsigma)
    ax.plot(geo, ones(length(geo)) * high, label="$nsigma sigma MAD cutoff")
    bigidx = means .> high 
    ax.scatter(geo[bigidx], means[bigidx], c=outlier_color)
    _remove_ticks(ax)
    _set_figsize(fig)  # Double the width here
    fig.savefig("outliers1.pdf", bbox_inches="tight", transparent=true, dpi=100)


    # 3. line  plot of insar after outliers vs GPS
    Blin = sum(InsarTimeseries.prepB(geolist17, intlist17), dims=2);
    before = (Blin \ unw_vals17) * p2c * 365 * 10  # Change to mm/year

    start_date=Date(2015,1,1)
    ref = nothing
    lw = 5
    fig, ax = plot_gps_los(station_name, before; ref=ref, lw=lw,
                           end_date=end_date, start_date=start_date,
                           ylim=(-3, 3), title="", bigfont=false, offset=false)

    # Now solve without outliers
    dts, los = get_gps_los(station_name, reference_station=ref, end_date=end_date, start_date=start_date)
    day_nums = _get_day_nums(dts)
    idxs = InsarTimeseries._good_idxs(geolist17[bigidx], intlist17)
    vals_after_outliers = unw_vals17[idxs]
    after = p2c * (Blin[idxs] \ vals_after_outliers)  # cm/day
    full_defo = after * (dts[end] - dts[1]).value
    # ax.plot(dts, -full_defo/2 .+ day_nums .* after, "c", lw=3, label="Insar")
    ax.plot(dts, day_nums .* after, "C0", lw=lw)
    _remove_ticks(ax)
    _set_figsize(fig, h, 2.3w)
    fig.savefig("outliers3.pdf", bbox_inches="tight", transparent=true, dpi=100)

    println("Before: $before mm/yr, $(before * (dts[end] - dts[1]).value /3650) total defo")
    println("After: $(after * 3650) mm/yr, $full_defo total defo")

    # 4. histogram before/after outlier
    fig, ax = plt.subplots()

    # ax = axes[ii, jj]
    # ax = axes[idx]
    vm = maximum(abs.(unw_vals17 * p2c))
    n1, _, _ = ax.hist(unw_vals17.* p2c, range=(-vm, vm), bins=bins, label="All data", color="red")
    n1, _, _ = ax.hist(vals_after_outliers .* p2c, range=(-vm, vm), bins=bins, label="Outliers removed", color="C0")
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers4.pdf", bbox_inches="tight", transparent=true, dpi=100)




    stds78 = Float64.(std.([h5read("gps_pixels_78.h5", name) for name in station_name_list78]))
    dists78 = [MapImages.latlon_to_dist(MapImages.station_latlon("TXKM"), MapImages.station_latlon(n)) for n in station_name_list78]
    df78 = DataFrame(dists=dists78, stds=stds78, names=station_name_list78)

    stds85 = Float64.(std.([h5read("gps_pixels_85.h5", name) for name in station_name_list85]))
    dists85 = [MapImages.latlon_to_dist(MapImages.station_latlon("TXKM"), MapImages.station_latlon(n)) for n in station_name_list85]
    df85 = DataFrame(dists=dists85, stds=stds85, names=station_name_list85)

    sm78 = lm(@formula(stds ~ dists), df78)
    a1, b1 = coef(sm78)
    sm85 = lm(@formula(stds ~ dists), df85)
    a2, b2 = coef(sm85)

    fig, ax = plt.subplots()
    ax.scatter(dists78, stds78, s=40)
    xs = 1:maximum(dists78)
    ax.plot(xs, a1 .+ b1 .* xs, lw=3)

    ax.scatter(dists85, stds85, s=40)
    xs = 1:maximum(dists85)
    ax.plot(xs, a2 .+ b2 .* xs, lw=3)

    ax.set_xlabel("km from TXKM")
    ax.set_ylabel("Std. Dev. of unw values")
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers5.pdf", bbox_inches="tight", transparent=true, dpi=100)

    return (geolist17, intlist17, unw_vals17)
end

function pnas_outlier_figure()
end


import StatsBase: countmap
import Images: label_components
_max_label(labels) = sort([tup for tup in countmap(vec(labels)) if tup[1] != 0], by=x->x[2], rev=true)[1][1]

function save_pnas_images(;years=[2016, 2017, 2018], 
                          fnames=["velocities_$(year)_current.h5" for year in years],
                          # fnames=["velocities_$(year)_linear_max800_sigma4.h5" for year in years],
                          outdir=".",
                          # zoomlats=(31.6, 30.9), zoomlons=(-103.8, -103.),  # Zoomed just to stripes
                          zoomlats=(31.9, 31.0), zoomlons=(-103.8, -102.85),  # Covers txkm/txmh
                          latbounds=(32.8, 30.7041),
                          dset="velos_shifted/1", cmap1="seismic_wide_y", cmap2="seismic_wide_y", # cmap2=rdylbl, 
                          vm1=10, vm2=10,
                          dem_mask_height=nothing)
    maxdays = maximum([_num_days(Sario.load_geolist_from_h5(f, dset)) for f in fnames])
    @show maxdays
    vm1 *= (maxdays / 3650)
    vm2 *= (maxdays / 3650)
    @show vm1, vm2 

    m1, m2 = nothing, nothing
    for f in fnames
        @show f
        days = _num_days(Sario.load_geolist_from_h5(f, dset))
        @show days

        m1 = MapImage(f, dset)
        m1[m1 .== 0] .= NaN    

        if !isnothing(dem_mask_height)
            dem = Sario.load("elevation_looked.dem")
            labels = label_components(dem .< 1200)
            maxl = _max_label(labels)
            m1[labels .!= maxl] .= NaN
        end

        m1 = m1[latbounds, (-110., -90.)] # no lon bounds really
        @show MapImages.grid_extent(m1)

        @show extremanan(m1)
        m1 .*= (days / 3650)  # Cumulative, in cm
        @show extremanan(m1)
        m2 = m1[zoomlats, zoomlons]
        out1 = joinpath(outdir, replace(f, ".h5"  => ".tif"))
        # out2 = replace(f, ".h5"  => "_inset.tif")
        out2 = joinpath(outdir, replace(f, ".h5"  => "_inset_gps.tif"))

        # Save as both png and tif
        save_img_geotiff(out1, m1, cmap=cmap1, vm=vm1)
        # save_img_geotiff(out2, m2, cmap=cmap2, vm=vm2)
        plt.imsave(replace(out1, ".tif" => ".png"), m1, cmap=cmap2, vmin=-vm2, vmax=vm2, dpi=150)
        plt.imsave(replace(out2, ".tif" => ".png"), m2, cmap=cmap2, vmin=-vm2, vmax=vm2, dpi=150)
    end
    # Now save one with a colorbar
    fig, axes = plt.subplots(1, 2)
    axim = axes[1].imshow(m1, cmap=cmap1, vmax=vm1, vmin=-vm1)
    fig.colorbar(axim, ax=axes[1])
    axim = axes[2].imshow(m2, cmap=cmap2, vmax=vm2, vmin=-vm2)
    fig.colorbar(axim, ax=axes[2])

    cbar_name = "colorbars.png"
    save_paper_figure(fig, cbar_name, false, 400)
end



function pnas_gps_error_table()
    p78 = "/data1/scott/pecos/path78-bbox2/igrams_looked/"
    p85 = "/data4/scott/path85/stitched/igrams_looked/"

   
    errors = Any[]
    push!(errors, all_stations)
    # errors85_4 = get_gps_error(p85*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)
    # errors78_4 = get_gps_error(p78*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)

    for f in ["velocities_2018_stackavg_max800.h5", "velocities_2018_current.h5", "velocities_2017_stackavg_max800.h5", "velocities_2017_current.h5"]
        fname = p85*f; @show fname
        push!(errors, round.(get_gps_error(fname, all_stations, dset="velos_shifted/1", verbose=false, window=7, shift=-0.), digits=2))
        fname = p78*f; @show fname
        push!(errors, round.(get_gps_error(fname, all_stations, dset="velos_shifted/1", verbose=false, window=7, shift=-0.), digits=2))
    end
    df = DataFrame(errors)
    show(stdout, MIME("text/latex"), df)
    @show extrema.(filter.(!isnan, errors[2:end]))
    @show rms.(errors[2:end])
    @show maximum.(filter.(!isnan, [abs.(a) for a in errors[2:end]]))
    return df
end

function write_gps_pixels(station_names, max_temp=800)
    geolist, intlist, valid_idxs = InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", max_temp)
    gpss = [get_stack_vals("unw_stack.h5", s, 3, "stack_flat_shifted", valid_idxs) for s in station_names]
    for (idx, s) in enumerate(station_names)
        h5write("gps_pixels.h5", s, gpss[idx])
    end
    Sario.save_intlist_to_h5("gps_pixels.h5", intlist)
    Sario.save_geolist_to_h5("gps_pixels.h5", geolist)
    Sario.save_dem_to_h5("gps_pixels.h5", Sario.load("dem.rsc"))
end
