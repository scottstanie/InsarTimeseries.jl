import CSV
using DataFramesMeta
import StatsBase: countmap
import Images: label_components
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
mdates = pyimport("matplotlib.dates")

include("./plotting.jl")
# Formatting functions:
_remove_ticks(ax) = ax.tick_params(
    axis = "both",
    which = "both",
    bottom = false,
    top = false,
    # labelbottom=false, labelleft=false,
    left = false,
    right = false,
)
_set_figsize(fig, h = 3, w = 3.5) = (fig.set_figheight(h), fig.set_figwidth(w))

function plot_rrc_oil_water()
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["font.family"] = "Helvetica"
    rcParams["font.size"] = 16
    rcParams["font.weight"] = "bold"

    data = readdf("summary_inj_prod.csv")
    oilyears = data[!, :production]
    wateryears = data[!, :injection]
    xs = data[!, :Year]

    fig, axes = plt.subplots(2, 1)
    # fig2, axes2 = plt.subplots(2, 1)
    bg_color = "b"
    active_color = "r"
    colors = fill(bg_color, length(xs))
    colors[end - 1] = active_color

    xticks = [2007, 2011, 2015, 2018]
    # yticks=[0, 0.5, 1, 1.5, 2, 2.5]

    titles = ["Oil Production", "Total Injections"]
    for (title, ax, values) in zip(titles, axes, [oilyears, wateryears])
        vals = values ./ 365 ./ 1e6
        # ax.bar(xs, vals, color=colors)
        ax.plot(xs, vals, "-o", lw = 5, ms = 10)
        ax.set_xticks(xticks)
        # ax.set_yticks(yticks)
        ax.set_ylabel("Million bbl/day")
        ax.grid()
        ax.set_title(title)
    end
    fig.set_figheight(9)
    fig.set_figwidth(4.5)
    fig.savefig("cisr-lines.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

end

function pnas_outlier_figure(
    station_name = "TXMC",
    outlier_color = "r",
    h = 3,
    w = 3.5,
    bins = 40,
)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # https://matplotlib.org/tutorials/introductory/customizing.html#a-sample-matplotlibrc-file
    # so obscure
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    rcParams["font.family"] = "Helvetica"
    rcParams["font.size"] = 16
    rcParams["font.weight"] = "bold"

    # a. histogram before/after outlier
    fig, ax = plt.subplots()

    unw_vals = h5read("gps_pixels_78.h5", station_name)
    geolist = Sario.load_geolist_from_h5("gps_pixels_78.h5")
    intlist = Sario.load_intlist_from_h5("gps_pixels_78.h5")
    end_date = Date(2018, 1, 1)
    geolist18 = geolist[geolist .> end_date]
    idxs = InsarTimeseries._good_idxs(geolist18, intlist)
    intlist17 = intlist[idxs]
    geolist17 = geolist[geolist .< end_date]
    unw_vals17 = unw_vals[idxs]
    unw_vals17_cm = p2c .* unw_vals17
    vm = maximum(abs.(unw_vals17_cm))
    n1, _, _ = ax.hist(
        unw_vals17_cm,
        range = (-vm, vm),
        bins = bins,
        label = "All data",
        color = "red",
    )
    MAT.matwrite("outliers_hist_vals.mat", Dict("unw_vals_$station_name" => unw_vals17_cm))
    # n1, _, _ = ax.hist(unw_vals17.* p2c, range=(-vm, vm), bins=bins, label="All data", color="red", density=true, alpha=0.5)
    ax.grid()

    means = InsarTimeseries.mean_abs_val(geolist17, intlist17, unw_vals17)
    nsigma = 4
    low, high = InsarTimeseries.two_way_cutoff(means, nsigma)
    bigidx = means .> high
    idxs = InsarTimeseries._good_idxs(geolist17[bigidx], intlist17)
    vals_after_outliers_cm = unw_vals17_cm[idxs]
    n1, _, _ = ax.hist(
        vals_after_outliers_cm,
        range = (-vm, vm),
        bins = bins,
        label = "Outliers removed",
        color = "C0",
    )
    # n1, _, _ = ax.hist(vals_after_outliers .* p2c, range=(-vm, vm), bins=bins, density=true, label="Outliers removed", color="C0", alpha=.5)
    ax.set_xlim((-15, 15))
    ax.set_xticks([-15, -7.5, 0, 7.5, 15])
    ax.set_yticks([0, 50, 100])
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers_a.pdf", bbox_inches = "tight", transparent = true, dpi = 100)


    # b) std dev of insar at GPS stations vs distance to reference
    dists = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in all_stations
    ]


    txmc_cm_78 = Float64.(p2c .* h5read("gps_pixels_78.h5", "TXMC"))
    @show median(txmc_cm_78)


    stds78_cm =
        Float64.(std.(p2c .* [h5read("gps_pixels_78.h5", name) for name in station_name_list78],))
    dists78 = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in station_name_list78
    ]

    stds85_cm =
        Float64.(std.(p2c .* [h5read("gps_pixels_85.h5", name) for name in station_name_list85],))
    dists85 = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in station_name_list85
    ]
    stations_all = vcat(station_name_list78, station_name_list85)
    dists_all = vcat(dists78, dists85)
    stds_all = vcat(stds78_cm, stds85_cm)
    dists_all = dists_all[stations_all .!= "TXKM"]
    stds_all = stds_all[stations_all .!= "TXKM"]
    stations_all = stations_all[stations_all .!= "TXKM"]


    df_dists = DataFrame(dists = dists_all, stds = stds_all, names = stations_all)

    sm = lm(@formula(stds ~ dists), df_dists)      # Linear model
    sm2 = lm(@formula(stds ~ sqrt(dists)), df_dists) # sqrt model
    # a, b = coef(sm)
    a, b = coef(sm2)
    @show a, b * 50

    fig, ax = plt.subplots()
    dists, stds = (df_dists[!, :dists], df_dists[!, :stds])
    ax.scatter(dists, stds, s = 40)
    xs = 1:maximum(dists)
    # ax.plot(xs, a .+ b .* xs, lw=3)  # Linear model
    ax.plot(xs, a .+ b .* sqrt.(xs), lw = 3)    # sqrt(dist) model
    ax.set_ylim((0, 5))
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_xticks([0, 100, 200])
    ax.set_ylabel("Noise Level [cm]")
    ax.set_xlabel("Distance [km]")
    ax.grid()
    # _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers_b.pdf", bbox_inches = "tight", transparent = true, dpi = 100)
    MAT.matwrite(
        "outliers_b.mat",
        Dict("a" => a, "b" => b, "xs" => collect(xs), "dists" => dists, "stds" => stds),
    )


    return sm, sm2
    # c) errors in 2, 3, 4 year solution at GPS vs distance, stacking only
    # d) errors in 2, 3, 4 year solution at GPS vs distance, after outlier removal
    df, errors = get_gps_errors(to_cm = true)
    df = @where(df, :station .!= "TXKM")

    years = collect(2016:2018)
    cols1 = [Symbol("p$(n)_stack_$(y)") for y in years, n in [78, 85]]
    cols2 = [Symbol("p$(n)_outliers_$(y)") for y in years, n in [78, 85]]
    colors = ["b", "g", "k"]
    fnames = ["outliers_c.pdf", "outliers_d.pdf"]

    maxy = 3
    for (idx, cols) in enumerate([cols1, cols2])
        fig, ax = plt.subplots()
        # errvals = abs.(vcat([df[!, c] for c in cols]...))
        for yidx = 1:3
            year = years[yidx]
            c1, c2 = cols[yidx, :]
            errvals = abs.(vcat(df[!, c1], df[!, c2]))
            dists = repeat(df[!, :dist_to_txkm], 2)
            dists = dists[.!isnan.(errvals)]
            errvals = errvals[.!isnan.(errvals)]
            @show maximum(errvals)

            ax.scatter(dists, errvals, s = 40, c = colors[yidx])
            ax.set_ylim((0, maxy))

            df1 = DataFrame(dist = dists, errval = errvals)
            sm = lm(@formula(errval ~ dist), df1)
            a, b = coef(sm)
            @show a, b * 50
            ax.plot(xs, a .+ b .* xs, colors[yidx], lw = 3)

            _remove_ticks(ax)
            _set_figsize(fig)
            fname = fnames[idx]
            fig.savefig(fname, bbox_inches = "tight", transparent = true, dpi = 100)
        end
    end
end


_max_label(labels) = sort(
    [tup for tup in countmap(vec(labels)) if tup[1] != 0],
    by = x->x[2],
    rev = true,
)[1][1]

function save_pnas_images(;
    years = [2016, 2017, 2018],
    fnames = ["velocities_$(year)_current.h5" for year in years],
    # fnames=["velocities_$(year)_linear_max800_sigma4.h5" for year in years],
    outdir = ".",
    # zoomlats=(31.6, 30.9), zoomlons=(-103.8, -103.),  # Zoomed just to stripes
    zoomlats = (31.9, 31.0),
    zoomlons = (-103.8, -102.85),  # Covers txkm/txmh
    latbounds = (32.8, 30.7041),
    dset = "velos_shifted/1",
    cmap1 = "seismic_wide_y",
    cmap2 = "seismic_wide_y", # cmap2=rdylbl, 
    vm1 = 10,
    vm2 = 10,
    dem_mask_height = nothing,
)
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

        m1 = m1[latbounds, (-110.0, -90.0)] # no lon bounds really
        @show MapImages.grid_extent(m1)

        @show extremanan(m1)
        m1 .*= (days / 3650)  # Cumulative, in cm
        @show extremanan(m1)
        m2 = m1[zoomlats, zoomlons]
        out1 = joinpath(outdir, replace(f, ".h5" => ".tif"))
        # out2 = replace(f, ".h5"  => "_inset.tif")
        out2 = joinpath(outdir, replace(f, ".h5" => "_inset_gps.tif"))

        # Save as both png and tif
        save_img_geotiff(out1, m1, cmap = cmap1, vm = vm1)
        # save_img_geotiff(out2, m2, cmap=cmap2, vm=vm2)
        plt.imsave(
            replace(out1, ".tif" => ".png"),
            m1,
            cmap = cmap2,
            vmin = -vm2,
            vmax = vm2,
            dpi = 150,
        )
        plt.imsave(
            replace(out2, ".tif" => ".png"),
            m2,
            cmap = cmap2,
            vmin = -vm2,
            vmax = vm2,
            dpi = 150,
        )
    end
    # Now save one with a colorbar
    fig, axes = plt.subplots(1, 2)
    axim = axes[1].imshow(m1, cmap = cmap1, vmax = vm1, vmin = -vm1)
    fig.colorbar(axim, ax = axes[1])
    axim = axes[2].imshow(m2, cmap = cmap2, vmax = vm2, vmin = -vm2)
    fig.colorbar(axim, ax = axes[2])

    cbar_name = "colorbars.png"
    save_paper_figure(fig, cbar_name, false, 400)
end


function _errs_to_cm(df)
    # mm/year, to cm for 2,3,4 years
    scales = vcat(4.1 .* ones(4), 3 .* ones(4), 2 .* ones(4)) ./ 10
    for (s, c) in zip(scales, names(df)[3:end])
        df[!, c] = df[!, c] .* s
    end
    return df
end


function get_gps_errors(; outfile = "gps_errors.csv", to_cm = false)
    if isfile(outfile)
        df = CSV.read(outfile)
        if to_cm
            df = _errs_to_cm(df)
        end
        return df, [c for c in eachcol(df)][3:end]
    end

    p78 = "/data1/scott/pecos/path78-bbox2/igrams_looked/"
    p85 = "/data4/scott/path85/stitched/igrams_looked/"


    dists = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in all_stations
    ]
    errors = Any[]
    push!(errors, all_stations)
    push!(errors, dists)
    # errors85_4 = get_gps_error(p85*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)
    # errors78_4 = get_gps_error(p78*"velocities_2018_current.h5", all_stations, dset="velos_shifted/1", verbose=true, window=7, shift=-0.)

    for f in [
        "velocities_2018_stackavg_max800.h5",
        "velocities_2018_current.h5",
        "velocities_2017_stackavg_max800.h5",
        "velocities_2017_current.h5",
        "velocities_2016_stackavg_max800.h5",
        "velocities_2016_current.h5",
    ]
        fname = p85 * f
        @show fname
        push!(
            errors,
            round.(get_gps_error(
                    fname,
                    all_stations,
                    dset = "velos_shifted/1",
                    verbose = false,
                    window = 7,
                    shift = -0.0,
                ),
                digits = 2,),
        )
        fname = p78 * f
        @show fname
        push!(
            errors,
            round.(get_gps_error(
                    fname,
                    all_stations,
                    dset = "velos_shifted/1",
                    verbose = false,
                    window = 7,
                    shift = -0.0,
                ),
                digits = 2,),
        )
    end
    columns = [
        :station,
        :dist_to_txkm,
        :p85_stack_2018,
        :p78_stack_2018,
        :p85_outliers_2018,
        :p78_outliers_2018,
        :p85_stack_2017,
        :p78_stack_2017,
        :p85_outliers_2017,
        :p78_outliers_2017,
        :p85_stack_2016,
        :p78_stack_2016,
        :p85_outliers_2016,
        :p78_outliers_2016,
    ]

    df = DataFrame(errors, columns)
    CSV.write("gps_errors.csv", df)
    return df, errors
end

function pnas_gps_error_table()
    df, errors = get_gps_errors()
    show(stdout, MIME("text/latex"), df)
    @show extrema.(filter.(!isnan, errors[3:end]))
    @show rms.(errors[3:end])
    @show maximum.(filter.(!isnan, [abs.(a) for a in errors[3:end]]))
    return df
end

function write_gps_pixels(station_names, max_temp = 800)
    geolist, intlist, valid_idxs =
        InsarTimeseries.load_geolist_intlist("unw_stack.h5", "geolist_ignore.txt", max_temp)
    gpss = [
        get_stack_vals("unw_stack.h5", s, 3, "stack_flat_shifted", valid_idxs)
        for s in station_names
    ]
    for (idx, s) in enumerate(station_names)
        h5write("gps_pixels.h5", s, gpss[idx])
    end
    Sario.save_intlist_to_h5("gps_pixels.h5", intlist)
    Sario.save_geolist_to_h5("gps_pixels.h5", geolist)
    Sario.save_dem_to_h5("gps_pixels.h5", Sario.load("dem.rsc"))
end

function pnas_outlier_figure_old(
    station_name = "TXMC",
    outlier_color = "r",
    h = 3,
    w = 3.5,
    bins = 40,
)

    # unw_vals = [h5read("gps_pixels.h5", name) for name in station_name_list78];
    unw_vals = h5read("gps_pixels_78.h5", station_name)
    geolist = Sario.load_geolist_from_h5("gps_pixels_78.h5")
    intlist = Sario.load_intlist_from_h5("gps_pixels_78.h5")

    end_date = Date(2018, 1, 1)
    geolist18 = geolist[geolist .> end_date]
    idxs = InsarTimeseries._good_idxs(geolist18, intlist)
    intlist17 = intlist[idxs]

    geolist17 = geolist[geolist .< end_date]
    unw_vals17 = unw_vals[idxs]
    geolist17, intlist17, unw_vals17 = (geolist, intlist, unw_vals)

    # 1. plot grouped by date
    fig, ax = plot_big_days(
        geolist17,
        intlist17,
        unw_vals17,
        legend = false,
        color = outlier_color,
    )
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers2.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

    # 2. cm vs baseline igram outlier plot
    # function plot_grouped_by_day(geo, int, vals, nsigma=0, min_spread=2)
    geo = geolist17
    means = InsarTimeseries.mean_abs_val(geo, intlist17, p2c * unw_vals17)
    println("median val: $(median(means))")
    fig, ax = plt.subplots()
    ax.scatter(geo, means, label = "means by date")
    nsigma = 4
    low, high = InsarTimeseries.two_way_cutoff(means, nsigma)
    ax.plot(geo, ones(length(geo)) * high, label = "$nsigma sigma MAD cutoff")
    bigidx = means .> high
    ax.scatter(geo[bigidx], means[bigidx], c = outlier_color)
    _remove_ticks(ax)
    _set_figsize(fig)  # Double the width here
    fig.savefig("outliers1.pdf", bbox_inches = "tight", transparent = true, dpi = 100)


    # 3. line  plot of insar after outliers vs GPS
    Blin = sum(InsarTimeseries.prepB(geolist17, intlist17), dims = 2)
    before = (Blin \ unw_vals17) * p2c * 365 * 10  # Change to mm/year

    start_date = Date(2015, 1, 1)
    ref = nothing
    lw = 5
    fig, ax = plot_gps_los(
        station_name,
        before;
        ref = ref,
        lw = lw,
        end_date = end_date,
        start_date = start_date,
        ylim = (-3, 3),
        title = "",
        bigfont = false,
        offset = false,
    )

    # Now solve without outliers
    dts, los = get_gps_los(
        station_name,
        reference_station = ref,
        end_date = end_date,
        start_date = start_date,
    )
    day_nums = _get_day_nums(dts)
    idxs = InsarTimeseries._good_idxs(geolist17[bigidx], intlist17)
    vals_after_outliers = unw_vals17[idxs]
    after = p2c * (Blin[idxs] \ vals_after_outliers)  # cm/day
    full_defo = after * (dts[end] - dts[1]).value
    # ax.plot(dts, -full_defo/2 .+ day_nums .* after, "c", lw=3, label="Insar")
    ax.plot(dts, day_nums .* after, "C0", lw = lw)
    _remove_ticks(ax)
    _set_figsize(fig, h, 2.3w)
    fig.savefig("outliers3.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

    println("Before: $before mm/yr, $(before * (dts[end] - dts[1]).value / 3650) total defo")
    println("After: $(after * 3650) mm/yr, $full_defo total defo")

    # 4. histogram before/after outlier
    fig, ax = plt.subplots()

    # ax = axes[ii, jj]
    # ax = axes[idx]
    vm = maximum(abs.(unw_vals17 * p2c))
    n1, _, _ = ax.hist(
        unw_vals17 .* p2c,
        range = (-vm, vm),
        bins = bins,
        label = "All data",
        color = "red",
    )
    n1, _, _ = ax.hist(
        vals_after_outliers .* p2c,
        range = (-vm, vm),
        bins = bins,
        label = "Outliers removed",
        color = "C0",
    )
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers4.pdf", bbox_inches = "tight", transparent = true, dpi = 100)




    stds78 =
        Float64.(std.([h5read("gps_pixels_78.h5", name) for name in station_name_list78]))
    dists78 = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in station_name_list78
    ]
    df78 = DataFrame(dists = dists78, stds = stds78, names = station_name_list78)

    stds85 =
        Float64.(std.([h5read("gps_pixels_85.h5", name) for name in station_name_list85]))
    dists85 = [
        MapImages.latlon_to_dist(
            MapImages.station_latlon("TXKM"),
            MapImages.station_latlon(n),
        ) for n in station_name_list85
    ]
    df85 = DataFrame(dists = dists85, stds = stds85, names = station_name_list85)

    sm78 = lm(@formula(stds ~ dists), df78)
    a1, b1 = coef(sm78)
    sm85 = lm(@formula(stds ~ dists), df85)
    a2, b2 = coef(sm85)

    fig, ax = plt.subplots()
    ax.scatter(dists78, stds78, s = 40)
    xs = 1:maximum(dists78)
    ax.plot(xs, a1 .+ b1 .* xs, lw = 3)

    ax.scatter(dists85, stds85, s = 40)
    xs = 1:maximum(dists85)
    ax.plot(xs, a2 .+ b2 .* xs, lw = 3)

    ax.set_xlabel("km from TXKM")
    ax.set_ylabel("Std. Dev. of unw values")
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("outliers5.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

    return (geolist17, intlist17, unw_vals17)
end

# Hunjoo figures
# 
function dict_to_vars(d::AbstractDict{String,Any})
    for key in keys(d)
        @eval $(Symbol(key)) = d[$(QuoteNode(key))]
    end
end
function get_wells()
    wolfidx = occursin.("WOLFCAMP", ff[!, prodcol])
    delidx = occursin.("DELAWARE", ff[!, prodcol])
end

function hunjoo_figures()
    rcParams["font.family"] = "Helvetica"
    rcParams["font.size"] = 16
    rcParams["font.weight"] = "medium"

    # variables = MAT.matread("/Users/scott/Documents/Learning/PNAS-git-overleaf/data/PNAS_Geomechanics3.mat")
    # variables = MAT.matread("/Users/scott/Desktop/figures/hunjoo/PNAS_Geomechanics7.mat")
    variables = MAT.matread("/Users/scott/Documents/Learning/GRL-git-overleaf-resubmit/hunjoo-revision-0526/GRL_rev1_Final_b0.01_Scott.mat")

    # dict_to_vars(variables)  # Fails cuz of some namespace thing
    # for key in keys(variables)
    #     @eval $(Symbol(key)) = variables[$(QuoteNode(key))]
    # end
    #  Area of Interest
    R_earth = 3959 * 5280 # (miles to feet) Earth Radius

    XD_min = -103.59
    XD_max = -103.44
    YD_min = 31.32
    YD_max = 31.47

    # #  Reservoir properties
    # RD = variables["Reservoir"][:, 1:2]   # RD = Reservoir(:,1),Reservoir(:,2)];
    # RR_R_D1 = (variables["Reservoir"][:, 3] / (5 * 3280.84) * 55).^2      # Reservoir radius to scatter plot (55^2 = 5 km = 5*3280.84 ft)
    # RR_R_D2 = (variables["Reservoir"][:, 3] / (1 * 3280.84) * 35).^2      # Reservoir radius to scatter plot (35^2 = 1 km = 3280.84 ft)

    # Cities
    Pecos = [-103.4932, 31.4229]             # (degrees)

    # Figure Scale
    ZD_min = -15.0
    ZD_max = 15.0
    c_min = -10.0
    c_max = 10.0
    FZ_min = -4.0
    FZ_max = 2.0

    demrsc = MapImages.from_grid(variables["XD"], variables["YD"])
    
    # Vertical Deformation
    cmap = "seismic_narrow_y"
    vm = 10

    # Pecos Area InSAR 2,3,1
    # TODO: DRAW TRANSECT FROM B B''
    fig, ax = plt.subplots()
    ax.imshow(variables["Z_Disp_VV18"], cmap=cmap, vmax=vm, vmin=-vm)
    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("hunjoo1.pdf", bbox_inches="tight", transparent=true, dpi=100)


    # Fault model: 2,3,2
    fig, ax = plt.subplots()
    # hold on
    # surf(XD_VV18,YD_VV18,Z_Disp_VV18,"EdgeColor","none")
    axim = ax.imshow(variables[ "Z_Disp_F3" ], cmap=cmap, vmax=vm, vmin=-vm, extent=(XD_min, XD_max,YD_min,YD_max))  #,"EdgeColor","none")   # InSAR data
    fig.colorbar(axim)
    ax.plot(variables[ "BB_XD" ][[1, end]], variables[ "BB_YD" ][[1, end]], lw=1.5)
    ax.plot(Pecos[1],Pecos[2], "ko");
    save_img_geotiff("hunjoo2.tif", MapImage(variables[ "Z_Disp_F3" ], demrsc), cmap = cmap, vm = vm)


    # reservoir model: 2,3,3
    fig, ax = plt.subplots()
    # hold on
    # surf(XD_VV18,YD_VV18,Z_Disp_VV18,"EdgeColor","none")
    cur_im = variables[ "Z_Disp_R3" ][:, :, 3]
    axim = ax.imshow(cur_im, cmap=cmap, vmax=vm, vmin=-vm, extent=(XD_min, XD_max,YD_min,YD_max))  #,"EdgeColor","none")   # InSAR data
    fig.colorbar(axim)
    ax.plot(variables[ "BB_XD" ][[1, end]], variables[ "BB_YD" ][[1, end]], lw=1.5)
    ax.plot(Pecos[1],Pecos[2], "ko")
    save_img_geotiff("hunjoo3.tif", MapImage(cur_im, demrsc), cmap = cmap, vm = vm)
    

    # Fault + reservoir model: 2,3,4
    fig, ax = plt.subplots()
    # surf(XD,YD,Z_Disp_T(:,:,3),"EdgeColor","none")
    im_mod = variables["Z_Disp_T3"][:, :, 3]
    axim = ax.imshow(
        im_mod,
        cmap = cmap,
        vmax = vm,
        vmin = -vm,
        extent = (XD_min, XD_max, YD_min, YD_max),
    )  # ,"EdgeColor","none")   # InSAR data
    # fig.colorbar(axim)
    ax.plot(variables["BB_XD"][[1, end]], variables["BB_YD"][[1, end]], lw = 1.5)
    ax.plot(Pecos[1], Pecos[2], "ko")

    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("hunjoo4.pdf", bbox_inches = "tight", transparent = true, dpi = 100)
    save_img_geotiff("hunjoo4.tif", MapImage(im_mod, demrsc), cmap = cmap, vm = vm)


    # Residual map: 2,3,5
    fig, ax = plt.subplots()
    # surf(XD,YD,Z_Disp_T(:,:,3),"EdgeColor","none")
    im_mod = variables["Res_Z_VV18"]
    axim = ax.imshow(
        im_mod,
        cmap = cmap,
        vmax = vm,
        vmin = -vm,
        extent = (XD_min, XD_max, YD_min, YD_max),
    )  # ,"EdgeColor","none")   # InSAR data
    # fig.colorbar(axim)
    ax.plot(variables["BB_XD"][[1, end]], variables["BB_YD"][[1, end]], lw = 1.5)
    ax.plot(Pecos[1], Pecos[2], "ko")

    _remove_ticks(ax)
    _set_figsize(fig)
    fig.savefig("hunjoo5.pdf", bbox_inches = "tight", transparent = true, dpi = 100)
    save_img_geotiff("hunjoo5.tif", MapImage(im_mod, demrsc), cmap = cmap, vm = vm)


    # Supplement: figure reservoir figure with circles
    fig, ax = plt.subplots()
    # surf(XD,YD,Z_Disp_T(:,:,3),"EdgeColor","none")
    im_mod = variables["UZ"]
    axim = ax.imshow(
        im_mod,
        cmap = cmap,
        vmax = vm,
        vmin = -vm,
        extent = (XD_min, XD_max, YD_min, YD_max),
    )
    _remove_ticks(ax)
    _set_figsize(fig)
    # fig.savefig("hunjoo-supp1.pdf", bbox_inches = "tight", transparent = true, dpi = 100)
    save_img_geotiff("hunjoo-supp1.tif", MapImage(im_mod, demrsc), cmap = cmap, vm = vm)


    
    # # B-B' cross-section
    # fig, ax1 = plt.subplots()
    # ax1.plot(BB_L', BB_Z_Disp_I18, "bs", label = "2018 InSAR")
    # ax1.plot(BB_L', BB_Z_Disp_R18, "r-.", lw = 1.5, label = "Reservoir")
    # ax1.plot(BB_L', BB_Z_Disp_F18, "r:", lw = 1.5, label = "Fault")
    # ax1.plot(BB_L', BB_Z_Disp_T18, "r-", lw = 1.5, label = "Res. + Fault")
    # ax1.set_xlim((0, maximum(BB_L)))
    # ax1.set_ylim((2 * FZ_min, 2 * FZ_max))
    # ax1.set_xlabel("B-B' (km)")
    # ax1.set_ylabel("Surface Deformation (cm)")

    # # yyaxis right
    # ax2 = ax1.twinx()
    # ax2.plot(BB_FC[:, 1, 1], BB_FC[:, 2, 1], "g-", lw = 2.0)
    # ax2.plot(BB_FC[:, 1, 2], BB_FC[:, 2, 2], "g-", lw = 2.0)
    # ax2.plot(BB_FC[:, 1, 3], BB_FC[:, 2, 3], "g-", lw = 2.0)
    # ax2.plot(BB_FC[:, 1, 4], BB_FC[:, 2, 4], "g-", lw = 2.0)
    # # text(BB_FC(2,1,1),BB_FC(2,2,1),"F1",lw=2.0)  #,"VerticalAlignment","top")
    # # text(BB_FC(2,1,2),BB_FC(2,2,2),"F2",lw=2.0)  #,"VerticalAlignment","top")
    # # text(BB_FC(2,1,3),BB_FC(2,2,3),"F3",lw=2.0)  #,"VerticalAlignment","top")
    # # text(BB_FC(2,1,4),BB_FC(2,2,4),"F4",lw=2.0)  #,"VerticalAlignment","top")
    # ax2.set_ylabel("Fault Depth (km)")
    # ax2.set_xlim((0, maximum(BB_L)))
    # ax2.set_ylim((0.5 * FZ_min, 0.5 * FZ_max))
    # # ax2.legend(["2018 InSAR","Reservoir Model","Fault Model","Res.+Fault Model","Location","northeast"])
    # # set(findobj(gcf,"type","axes"),"FontSize",12,"TickDir","out")
    # #
    # xticks = [0, 7, 14]
    # y1ticks = [-8, -4, 0, 4]
    # y2ticks = [-2, -1, 0, 1]

    # ax1.set_xticks(xticks)
    # ax1.set_yticks(y1ticks)
    # ax2.set_yticks(y2ticks)
    # ax1.grid()
    # fig.legend(frameon = false, fontsize = "small", ncol = 2)

    # _set_figsize(fig, 1.3 .* [3, 3.5]...)
    # # _set_figsize(fig)
    # fig.savefig("hunjoo4.pdf", bbox_inches = "tight", transparent = true, dpi = 100)
end


# Functions for weighting by covariance

import ImageFiltering: imfilter, Kernel
function get_aps(stack_in)
    return imfilter(stack_in .- imfilter(stack_in, (Kernel.gaussian((0, 0, 1)), Kernel.Laplacian((false, false, true)))),  Kernel.gaussian((2, 2, 0)))
end

function try_weight()
    test_sub = get_aps(defo_stack_unreg);
    var_by_day = dropdims(var((defo_stack_unreg .- test_sub); dims = (1, 2)), dims = (1, 2));
    var_by_day2 = abs.(test_sub[400, 400, :])
    Qddfull2 = build_Q(abs.(var_by_day2), geolist, intlist);
    Qinv2 = inv(Qddfull2);
    W = inv(A' * Qinv * A) * A' * Qinv2;
end

function build_Q(var_dates, geolist, intlist, diag_val = 0.5)
    N = length(intlist);
    Q = zeros((N, N));
    for (jdx, ig2) in enumerate(intlist)
        for (idx, ig1) in enumerate(intlist)
            if jdx > idx
            # if Q[jdx, idx] != 0 
                Q[idx, jdx] = Q[jdx, idx]
                continue  # dont need to redo, symmetric
            end
            # @show idx, jdx
            if idx == jdx
                v1 = var_dates[findfirst(isequal(ig1[1]), geolist)]
                v2 = var_dates[findfirst(isequal(ig1[2]), geolist)]
                Q[idx, jdx] = v1 + v2
            else
                cm = countmap([ig1..., ig2...])
                if length(Set([ig1..., ig2...])) == 4
                    Q[idx, jdx] = 0
                    continue
                end
                repeated_date = sort(collect(countmap([ig1..., ig2...])), by = x->x[2], rev = true)[1][1]
                # @show repeated_date
                # @show findfirst(isequal(repeated_date), geolist)
                sigma_repeat = var_dates[findfirst(isequal(repeated_date), geolist)]
                # if jdx in (1, 2) && idx == 28
                    # @show idx, jdx, repeated_date, sigma_repeat
                if (repeated_date == ig1[1] && repeated_date == ig2[1]) || (repeated_date == ig1[2] && repeated_date == ig2[2])
                    Q[idx, jdx] = sigma_repeat
                else
                    Q[idx, jdx] = -sigma_repeat
                end
            end
        end
    end
    return Q + diag_val * I
end


_intlsubset(geolist, intlist, N) = intlist[[(intl[1] in geolist[1:N] && intl[2] in geolist[1:N]) for intl in intlist]]
function _intlsubset(geolist, intlist, vallist::AbstractArray, N)
    idxs = [(intl[1] in geolist[1:N] && intl[2] in geolist[1:N]) for intl in intlist]
    intlist[idxs], vallist[idxs]
end

function Qsubset(var_by_day, geolist, intlist, N)
    g = geolist[1:N]
    v = var_by_day[1:N]
    igl = _intlsubset(geolist, intlist, N)
    return build_Q(v, g, igl)
end

function build_Q_diag(var_dates, geolist, intlist)
    N = length(intlist);
    Q1 = zeros((N,));
    for (idx, ig) in enumerate(intlist)
        v1 = var_dates[findfirst(isequal(ig[1]), geolist)]
        v2 = var_dates[findfirst(isequal(ig[2]), geolist)]
        Q1[idx] = v1 + v2
    end
    return diagm(Q1)
end

import Combinatorics
function demo_weighted(bad_var = 100)
    geolist1 = Date(2020, 1, 1) .+ Dates.Day(12) .* (0:3)
    # 4-element Array{Date,1}:
    #  2020-01-01
    #  2020-01-13
    #  2020-01-25
    #  2020-02-06
    intlist1 = [tuple(p...) for p in Combinatorics.combinations(geolist1, 2)]
    # 6-element Array{Tuple{Date,Date},1}:
    #  (2020-01-01, 2020-01-13)
    #  (2020-01-01, 2020-01-25)
    #  (2020-01-01, 2020-02-06)
    #  (2020-01-13, 2020-01-25)
    #  (2020-01-13, 2020-02-06)
    #  (2020-01-25, 2020-02-06)

    # Real phase is 0, 1, 2, 3
    phi_true = [0, 1, 2, 3] 
    A = InsarTimeseries.build_A_matrix(geolist1, intlist1)
    @show meas_nonoise = A * phi_true[2:end]
    @show A \ meas_nonoise

    # But with noisy second day, looks like it uplifted 10
    phi_noise = [0, 10, 2, 3] 
    @show meas_noise = A * phi_noise[2:end]
    println("No weighting:")
    @show A \ meas_noise

    variances = [1, bad_var, 1, 1]
    a, b, c, d = variances
    # Q = [a + b a -b; a a + c c; -b c b + c] + 0.1 * I
    Q = build_Q(variances, geolist1, intlist1, .1)
    print("Weight LS:")
    @show inv(A' * inv(Q) * A) * A' * inv(Q) * meas_noise

    println("or just diagonal?")
    Q2 = build_Q_diag(variances, geolist1, intlist1)
    println(Q2)
    print("Weight LS:")
    @show inv(A' * inv(Q2) * A) * A' * inv(Q2) * meas_noise

end


function plot_l1_vs_stack(;offset=true, alpha=300, h=3, w=4.5, yy=2018, 
                          unwfile="unw_stack_shiftedonly.h5", station="TXOE", s=station,
                          days_smooth=10)
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    # https://matplotlib.org/tutorials/introductory/customizing.html#a-sample-matplotlibrc-file
    # so obscure
    rcParams["pdf.fonttype"] = 42
    rcParams["ps.fonttype"] = 42
    rcParams["font.family"] = "Helvetica"
    rcParams["font.size"] = 16
    rcParams["font.weight"] = "bold"
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    years_fmt = mdates.DateFormatter("%Y")
    # yrs = (Date(2015), Date(2016), Date(2017), Date(2018))

    demrsc = Sario.load("dem.rsc")
    input_dset = "stack_flat_shifted"
    geolist, intlist, igram_idxs = load_geolist_intlist(unwfile, "geolist_ignore.txt", 800, max_date = Date(yy,1,1))
    Blin = sum(InsarTimeseries.prepB(geolist, intlist), dims = 2)
    ms = 4
    timediffs = InsarTimeseries.day_diffs(geolist)
    unw_vals = get_stack_vals(unwfile, gps.station_rowcol(s, Dict(demrsc))..., 1, input_dset, igram_idxs)

    # Timeseries: unregularized, with all outliers (noisiest)
    B = InsarTimeseries.prepB(geolist, intlist, false, 0)
    unregged = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(B \ unw_vals, timediffs)

    # Timeseries: regularized, but with all outliers
    Ba = InsarTimeseries.prepB(geolist, intlist, false, alpha)
    unw_vals_a = InsarTimeseries.augment_zeros(Ba, unw_vals)
    regged = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(Ba \ unw_vals_a, timediffs)

    # 
    # Plot timeseries with-outlier cases: 
    fig, ax = plot_gps_los(s; end_date = Date(yy, 1, 1), offset=offset, 
                           days_smooth = days_smooth,
                           gps_color=matlab_colors()[4], ms=ms)

    ax.plot(geolist, unregged, "-x", lw=3, label="Regularized")
    ax.plot(geolist, regged, "-x", lw=3, label="Unregularized")
    ax.format_xdata = years_fmt
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    # ax.tick_params(axis="x", labelrotation=45)
    y0 = ceil(maximum(abs.(unregged)))
    ax.set_ylim((-y0, y0))
    ax.set_yticks(-y0:2:y0)
    _set_figsize(fig, h, w)
    fig.legend()
    fig.savefig("compare_$(yy)_timeseries.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

    # No outlier removal linear cases
    stack, l2, l1 = prunesolve(geolist, intlist, unw_vals, Blin, 1000, shrink = false)
    println("No outlier, linear: $(stack), $(l1)")
    println("Difference: $(abs(stack - l1))")


    # Plot linear with-outlier cases
    fig, ax = plot_gps_los(s; end_date = Date(yy, 1, 1), 
                             insar_mm_list=[l1, l2], offset=offset, 
                             labels=["L1 linear", "L2 linear"], 
                             gps_color=matlab_colors()[3],
                             insar_colors=matlab_colors()[1:2],
                             days_smooth = days_smooth,
                             ms=ms)

    # ax.plot(geolist, regged, "-x", lw=3, label="reg")
    ax.plot(geolist, unregged, "-x", lw=3, c=matlab_colors()[4], label="Unregularized")
    ax.plot(geolist, regged, "-x", lw=3, c=matlab_colors()[5], label="Regularized")
    ax.format_xdata = years_fmt
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    # ax.tick_params(axis="x", labelrotation=45)
    # y0 = ceil(maximum(abs.(regged)))
    y0 = 6
    ax.set_ylim((-y0, y0))
    ax.set_yticks(-y0:2:y0)
    _set_figsize(fig, h, w)
    fig.legend()
    fig.savefig("compare_$(yy)_linear.pdf", bbox_inches = "tight", transparent = true, dpi = 100)

    ###### Outlier remove cases ################
    # timeseries: regularized with outliers removed
    geo_clean, intlist_clean, unw_clean = remove_outliers(geolist, intlist, unw_vals, mean_sigma_cutoff = 4)
    B2 = InsarTimeseries.prepB(geo_clean, intlist_clean, false, 0)
    td_clean = InsarTimeseries.day_diffs(geo_clean)
    unregged2 = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(B2 \ unw_clean, td_clean)
    Ba2 = InsarTimeseries.prepB(geo_clean, intlist_clean, false, alpha)
    unw_vals_a2 = InsarTimeseries.augment_zeros(Ba2, unw_clean)
    regged2 = InsarTimeseries.PHASE_TO_CM .* InsarTimeseries.integrate_velocities(Ba2 \ unw_vals_a2, td_clean)


    # linear solves with outlier removal
    stack2, l22, l12 = prunesolve(geolist, intlist, unw_vals, Blin, 4, shrink = false)

    # PLOT:
    fig, ax = plot_gps_los(s; end_date = Date(yy, 1, 1),
                             insar_mm_list=[l12, l22], offset=offset, 
                             labels=["L1 linear", "L2 linear"], 
                             gps_color=matlab_colors()[3],
                             insar_colors=matlab_colors()[1:2],
                             days_smooth = days_smooth,
                             ms=ms)
    # ax.plot(geo_clean, regged2, "-x", lw=3, label="reg")
    ax.plot(geo_clean, unregged2, "-x", c=matlab_colors()[4], lw=3, label="Unregularized")
    ax.plot(geo_clean, regged2, "-x", lw=3, c=matlab_colors()[5], label="Regularized")
    # y0 = ceil(maximum(abs.(regged2)))
    y0 = 4
    ax.set_ylim((-y0, y0))
    ax.set_yticks(-y0:2:y0)
    ax.format_xdata = years_fmt
    ax.xaxis.set_major_locator(years)
    ax.xaxis.set_major_formatter(years_fmt)
    _set_figsize(fig, h, w)
    fig.legend()
    fig.savefig("compare_$(yy)_removed.pdf", bbox_inches = "tight", transparent = true, dpi = 100)


end