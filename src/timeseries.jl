using Dates
using Statistics
using LinearAlgebra
import TOML

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / ( 4 * π )
const P2MM = 365 * 10 * PHASE_TO_CM   # mm / year

const DateOrNone = Union{Date,Nothing}

_stack_size_mb(filename, dset) =
    prod(size(filename, dset)) * sizeof(eltype(filename, dset)) / 1e6
_strnothing(aa) = isnothing(aa) ? "nothing" : aa   # Since we can't string(nothing)


function run_inversion(config_file::Dict{AbstractString,Any})
    conf_dict = TOML.parsefile(config_file)
    # Convert to symbols so it works to pass as kwargs
    symbol_dict = Dict(Symbol(k) => v for (k, v) in conf_dict)
    outfile, outdset = run_inversion(symbol_dict...)
    h5writeattr(outfile, outdset, conf_dict)
end

# overwrites in the dset fail now... is that fine? force user to delete?
function run_inversion(;
    unw_stack_file::String=UNW_FILENAME,
    input_dset::String=STACK_FLAT_SHIFTED_DSET,
    outfile::String="",
    # outgroup::String="velos",
    stack_average::Bool=false,
    constant_velocity::Bool=true,
    ignore_geo_file::String="",
    max_temporal_baseline::Int=500,
    prune_outliers=true,
    sigma=4,
    gap=1,
    min_date::DateOrNone=nothing,
    max_date::DateOrNone=nothing,
    alpha::Real=0.0,
    L1::Bool=false,
    reference_station=nothing, # Note: these only work for stackavg...
    ref_row=nothing,
    ref_col=nothing,
    use_distributed=true,
)
    if isempty(outfile)
        outfile = _default_outfile()
    end

    # averaging or linear means output will is 3D array (not just map of velocities)
    is_3d = !(stack_average || constant_velocity)
    outdset = is_3d ? "stack" : "velos"

    # Check if this file/dset already exists
    isfile(outfile) &&
    outdset in names(outfile) &&
    error("$outdset already exists in $outfile")

    # the valid igram indices is out of all possible layers in the stack
    geolist, intlist, valid_igram_indices = load_geolist_intlist(
        unw_stack_file,
        ignore_geo_file,
        max_temporal_baseline,
        min_date=min_date,
        max_date=max_date,
    )


    println("$outdset geolist range: : $(extrema(geolist))")
    println("Using $input_dset as input")
    # Now: can we load the input stack into memory? or do we need distributed?
    stack_file_size = _stack_size_mb(unw_stack_file, input_dset)
    can_fit_mem = (stack_file_size * 8) < getmemavail()  # Rough padding for total memory check

    # Note: as of version 1.3.0rc2, my threaded version on preloaded matrix is slower than
    # distributed reading/writing... not sure why but oh well
    # TODO: figure if it's worth having non-process version
    println("Stack size: $stack_file_size, avail RAM: $(getmemavail()), fitting in ram: $can_fit_mem")
    println("Using Distributed to solve: $use_distributed")

    if stack_average
        # # Dont need shift for avg
        # input_dset = STACK_FLAT_DSET
        println("Averaging stack for solution")
        velo_file_out = run_stackavg(
            unw_stack_file,
            input_dset,
            outfile,
            outdset,
            geolist,
            intlist;
            reference_station=reference_station,
            ref_row=ref_row,
            ref_col=ref_col,
        )
    else
        println("Running SBAS for solution")

        velo_file_out = run_sbas(
            unw_stack_file,
            input_dset,
            outfile,
            outdset,
            geolist,
            intlist,
            valid_igram_indices,
            constant_velocity,
            alpha,
            L1,
            prune_outliers,
            sigma,
        )
        # if use_distributed
        # else
        #     # If we want to read the whole stack in at once:
        #     @time unw_stack = h5read(unw_stack_file, input_dset)[:, :, valid_igram_indices]
        #     velo_file_out = run_sbas(unw_stack, outfile, outdset, geolist, intlist, 
        #                              constant_velocity, alpha, L1, prune)
        # end

    end

    ####################

    dem_rsc = Sario.load_dem_from_h5(unw_stack_file)
    Sario.save_dem_to_h5(outfile, dem_rsc, overwrite=true)
    Sario.save_geolist_to_h5(outfile, outdset, geolist, overwrite=true)
    save_reference(outfile, unw_stack_file, outdset, input_dset)

    return outfile, outdset
end


function load_geolist_intlist(
    unw_stack_file,
    ignore_geo_file,
    max_temporal_baseline;
    min_date::DateOrNone=nothing,
    max_date::DateOrNone=nothing,
)
    @show unw_stack_file
    geolist = Sario.load_geolist_from_h5(unw_stack_file)
    intlist = Sario.load_intlist_from_h5(unw_stack_file)

    # If we are ignoreing some indices, remove them from for all pixels
    geo_idxs, igram_idxs = find_valid_indices(
        geolist,
        intlist,
        min_date,
        max_date,
        ignore_geo_file,
        max_temporal_baseline,
    )
    return geolist[geo_idxs], intlist[igram_idxs], igram_idxs
end

"""Cut down the full list of interferograms and geo dates

- Cuts date ranges (e.g. piecewise linear solution) with  `min_date` and `max_date`
- Can ignore whole dates by reading `ignore_geo_file`
- Prunes long time baseline interferograms with `max_temporal_baseline`
"""
function find_valid_indices(
    geo_date_list::Array{Date,1},
    igram_date_list::Array{Igram,1},
    min_date::DateOrNone=nothing,
    max_date=DateOrNone = nothing,
    ignore_geo_file::String="",
    max_temporal_baseline::Int=500,
)

    ig1 = length(igram_date_list)  # For logging purposes, what do we start with
    if isempty(ignore_geo_file)
        println("Not ignoring any .geo dates")
        ignore_geos = []
    else
        ignore_geos = sort(Sario.find_geos(filename=ignore_geo_file, parse=true))
        println("Ignoring the following .geo dates:")
        println(ignore_geos)
    end

    # TODO: do I want to be able to pass an array of dates to ignore?

    # First filter by remove igrams with either date in `ignore_geo_file`
    valid_geos = [g for g in geo_date_list if !(g in ignore_geos)]
    valid_igrams = [ig for ig in igram_date_list if !(ig in ignore_geos)]
    println("Ignoring $(ig1 - length(valid_igrams)) igrams listed in $ignore_geo_file")

    # Remove geos and igrams outside of min/max range
    if !isnothing(min_date)
        println("Keeping data after min_date: $min_date")
        valid_geos = [g for g in valid_geos if g > min_date]
        valid_igrams = [ig for ig in valid_igrams if (ig[1] > min_date && ig[2] > min_date)]
    end
    if !isnothing(max_date)
        println("Keeping data only before max_date: $max_date")
        valid_geos = [g for g in valid_geos if g < max_date]
        valid_igrams = [ig for ig in valid_igrams if (ig[1] < max_date && ig[2] < max_date)]
    end

    # This is just for logging purposes:
    too_long_igrams =
        filter(ig -> temporal_baseline(ig) > max_temporal_baseline, valid_igrams)
    println("Ignoring $(length(too_long_igrams)) igrams with longer baseline than $max_temporal_baseline days")

    ### Remove long time baseline igrams ###
    valid_igrams =
        filter(ig -> temporal_baseline(ig) <= max_temporal_baseline, valid_igrams)


    ### Collect remaining geo dates and igrams
    valid_geo_indices = indexin(valid_geos, geo_date_list)
    valid_igram_indices = indexin(valid_igrams, igram_date_list)

    println("Ignoring $(ig1 - length(valid_igrams)) igrams total")

    return valid_geo_indices, valid_igram_indices
end

"""Get min,max pairs to split up interval of dates
gap=1 means form adjacent pairs
gap=2 skips a date in between, makes overlapping pairs to get larger, redundant estimates
"""
function _get_pairs(dates, gap=1)
    isempty(dates) && return (nothing, nothing)
    length(dates) < gap && error("Need at least $gap dates if gap=$gap")
    dates_pad = [nothing; dates; nothing]
    return collect(zip(dates_pad[1:end - gap], dates_pad[(gap + 1):end]))
end



function save_reference(
    h5file,
    unw_stack_file,
    dset_name,
    stack_flat_shifted_dset,
    overwrite=true,
)
    # Delete if exists
    if overwrite
        h5open(h5file) do f
            if exists(attrs(f[dset_name]), REFERENCE_ATTR)
                a_delete(f[dset_name], REFERENCE_ATTR)
            end
            if exists(attrs(f[dset_name]), REFERENCE_STATION_ATTR)
                a_delete(f[dset_name], REFERENCE_STATION_ATTR)
            end
        end
    end
    # Read ref. unfo from unw file, save what was used to deformation result file
    reference = get(h5readattr(unw_stack_file, stack_flat_shifted_dset), REFERENCE_ATTR, "")
    h5writeattr(h5file, dset_name, Dict(REFERENCE_ATTR => reference))
    reference_station =
        get(h5readattr(unw_stack_file, stack_flat_shifted_dset), REFERENCE_STATION_ATTR, "")
    h5writeattr(h5file, dset_name, Dict(REFERENCE_STATION_ATTR => reference_station))
end
