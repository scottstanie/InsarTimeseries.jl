using Dates
using Statistics
using LinearAlgebra
import TOML

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / (-4 * π )
const P2MM = 365 * 10 * PHASE_TO_CM   # mm / year

const DateOrNone = Union{Date, Nothing}

_stack_size_mb(filename, dset) = prod(size(filename, dset)) * sizeof(eltype(filename, dset)) / 1e6 

"""Runs SBAS inversion on all unwrapped igrams

"""
function run_inversion(config_file::Dict{AbstractString, Any})
    conf_dict = TOML.parsefile(config_file)
    # Convert to symbols so it works to pass as kwargs
    symbol_dict = Dict(Symbol(k) => v for (k, v) in conf_dict)
    return run_inversion(symbol_dict...)
end 

function run_inversion(; unw_stack_file::String="",
                       input_dset::String=STACK_FLAT_SHIFTED_DSET,
                       outfile::String="", 
                       outdset::String="velos",
                       stack_average::Bool=false, 
                       constant_velocity::Bool=true, 
                       ignore_geo_file::String="", 
                       max_temporal_baseline::Int=500,
                       split_dates::AbstractArray=[],
                       split_count=1,
                       min_date::DateOrNone=nothing,
                       max_date::DateOrNone=nothing,
                       alpha::Float32=0.0f0,
                       L1::Bool=false,  
                       use_distributed=true,
                       kwargs...)
    if isempty(outfile)
        outfile = _default_outfile()
    end

    # Now for each split date, run this function on a section
    if !isempty(split_dates)
        for (d1, d2) in _get_pairs(split_dates)
            return run_inversion(split_dates=[], unw_stack_file=unw_stack_file, 
                                 input_dset=input_dset, outfile=outfile, outdset=outdset, 
                                 stack_average=stack_average, constant_velocity=constant_velocity, 
                                 ignore_geo_file=ignore_geo_file, 
                                 max_temporal_baseline=max_temporal_baseline, alpha=alpha, L1=L1, 
                                 use_distributed=use_distributed, kwargs...)
        end
    end
    
    # the valid igram indices is out of all possible layers in the stack
    geolist, intlist, valid_igram_indices = load_geolist_intlist(unw_stack_file, ignore_geo_file, 
                                                                 max_temporal_baseline,
                                                                 min_date=min_date,
                                                                 max_date=max_date)

    # # Dont need shift for avg
    # input_dset = stack_average ? STACK_FLAT_DSET : STACK_FLAT_SHIFTED_DSET
    # Now: can we load the input stack into memory? or do we need distributed?
    stack_file_size = _stack_size_mb(unw_stack_file, input_dset)
    can_fit_mem = (stack_file_size * 8) < getmemavail()  # Rough padding for total memory check

    # Note: as of version 1.3.0rc2, my threaded version on preloaded matrix is slower than
    # distributed reading/writing... not sure why but oh well
    println("Stack size: $stack_file_size, avail RAM: $(getmemavail()), fitting in ram: $can_fit_mem")
    println("Using Distributed to solve: $use_distributed")

    # outdset = stack_average ? "stack" : "velos"  # TODO: get this back to hwere not only velos
    if stack_average
        is_hdf5 = false
        println("Averaging stack for solution")
        vstack = run_stackavg(unw_stack_file, input_dset, geolist, intlist)
        is_3d = true  # TODO: stack avg should really be 2d velo
    else
        is_hdf5 = true
        println("Performing SBAS solution")
        println("Reading unw stack")
        # @time unw_stack = load_hdf5_stack(unw_stack_file, STACK_FLAT_SHIFTED_DSET)
        # @time unw_stack = load_hdf5_stack(unw_stack_file, STACK_FLAT_SHIFTED_DSET, valid_igram_indices)
        # TODO: also do this for the masks

        # vstack = run_sbas(unw_stack, geolist, intlist, constant_velocity, alpha)
        # outdset = "stack"
        # h5open(unw_stack_file) do unw_file
        #     unw_stack = unw_file[input_dset]
        #     vstack = run_sbas(unw_stack, geolist, intlist, valid_igram_indices,
        #                       constant_velocity, alpha, L1)
        # end
        if use_distributed
            velo_file_out = run_sbas(unw_stack_file, input_dset, outfile, outdset,
                                     geolist, intlist, valid_igram_indices, 
                                     constant_velocity, alpha, L1)
        else
            # If we want to read the whole stack in at once:
            @time unw_stack = h5read(unw_stack_file, input_dset)[:, :, valid_igram_indices]
            velo_file_out = run_sbas(unw_stack, outfile, outdset, geolist, intlist, 
                                     constant_velocity, alpha, L1)
        end

        is_3d = false
    end
    ####################33

    if is_3d
        timediffs = day_diffs(geolist)
        println("Integrating velocities to phases")
        @time phi_arr = integrate_velocities(vstack, timediffs)
        # Multiply by wavelength ratio to go from phase to cm
        deformation = PHASE_TO_CM .* phi_arr
    end

    dem_rsc = sario.load_dem_from_h5(unw_stack_file) 

    if is_3d
        println("Saving deformation to $outfile")
        @time save_deformation(outfile, deformation, dem_rsc, unw_stack_file=unw_stack_file, do_permute=!is_hdf5)
    else
        sario.save_dem_to_h5(outfile, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=true)
    end

    save_geolist_to_h5(outfile, geolist, overwrite=true)
    save_reference(outfile, unw_stack_file, outdset, input_dset)

    return outfile, outdset
end


function load_geolist_intlist(unw_stack_file, ignore_geo_file, max_temporal_baseline;
                             min_date::DateOrNone=nothing, max_date::DateOrNone=nothing)
    geolist = load_geolist_from_h5(unw_stack_file)
    intlist = load_intlist_from_h5(unw_stack_file)

    # If we are ignoreing some indices, remove them from for all pixels
    geo_idxs, igram_idxs = find_valid_indices(geolist, intlist, min_date, max_date, 
                                              ignore_geo_file, max_temporal_baseline)
    return geolist[geo_idxs], intlist[igram_idxs], igram_idxs
end

"""Cut down the full list of interferograms and geo dates

- Cuts date ranges (e.g. piecewise linear solution) with  `min_date` and `max_date`
- Can ignore whole dates by reading `ignore_geo_file`
- Prunes long time baseline interferograms with `max_temporal_baseline`
"""
function find_valid_indices(geo_date_list::Array{Date, 1}, igram_date_list::Array{Igram, 1}, 
                            min_date::DateOrNone=nothing, max_date=DateOrNone=nothing,
                            ignore_geo_file::String="", max_temporal_baseline::Int=500)

    if isempty(ignore_geo_file)
        println("Not ignoring any .geo dates")
    else
        ignore_geos = sort(sario.find_geos(filename=ignore_geo_file, parse=true))
        println("Ignoring the following .geo dates:")
        println(ignore_geos)
    end

    # TODO: do I want to be able to pass an array of dates to ignore?

    # First filter by remove igrams with either date in `ignore_geo_file`
    valid_geos = [g for g in geo_date_list if !(g in ignore_geos)]
    valid_igrams = [ig for ig in igram_date_list if !(ig in ignore_geos)]

    # Remove geos and igrams outside of min/max range
    if !isnothing(min_date)
        println("Removing data before $min_date")
        valid_geos = [g for g in valid_geos if  g < min_date]
        valid_igrams = [ig for ig in valid_igrams if (ig[1] < min_date || ig[2] < min_date)]
    end
    if !isnothing(max_date)
        println("Removing data after $max_date")
        valid_geos = [g for g in valid_geos if g > max_date ]
        valid_igrams = [ig for ig in valid_igrams if (ig[1] > max_date || ig[2] > max_date)]
    end

    ### Remove long time baseline igrams ###
    valid_igrams = filter(ig -> temporal_baseline(ig) <= max_temporal_baseline, valid_igrams)

    # This is just for logging purposes:
    too_long_igrams = filter(ig -> temporal_baseline(ig) > max_temporal_baseline, valid_igrams)
    println("Ignoring $(length(too_long_igrams)) igrams with longer baseline than $max_temporal_baseline days")


    ### Collect remaining geo dates and igrams
    valid_geo_indices = indexin(valid_geos, geo_date_list)
    valid_igram_indices = indexin(valid_igrams, igram_date_list)

    println("Ignoring $(length(igram_date_list) - length(valid_igrams)) igrams total")

    return valid_geo_indices, valid_igram_indices
end

"""Finds the number of days between successive .geo files"""
function day_diffs(geolist::Array{Date, 1})
    [difference.value for difference in diff(geolist)]
end


"""Takes velocity solution output and finds phases

Arguments:
    velocities come from invert_sbas
    timediffs are the days between each SAR acquisitions
        length will be 1 less than num SAR acquisitions
"""
function integrate_velocities(vstack::Array{<:AbstractFloat, 3}, timediffs::Array{Int, 1})
    nrows, ncols, _ = size(vstack)
    num_geos = length(timediffs) + 1
    phi_stack = zeros(nrows, ncols, num_geos)
    phi_diffs = zeros(num_geos - 1)

    phi_arr = zeros(num_geos)  # Buffer to hold each result
    for j = 1:ncols
        for i = 1:nrows
            varr = view(vstack, i, j, :)
            phi_stack[i, j, :] .= integrate1D!(phi_diffs, phi_arr, varr, timediffs)
        end
    end
    return phi_stack
end

function integrate1D!(phi_diffs, phi_arr, velocities::AbstractArray{<:AbstractFloat, 1}, timediffs::Array{Int, 1})
    # multiply each column of vel array: each col is a separate solution
    phi_diffs .= velocities .* timediffs

    # Now the final phase results are the cumulative sum of delta phis
    # This is equivalent to replacing missing with 0s (like np.ma.cumsum does)
    phi_arr[2:end] .= cumsum(coalesce.(phi_diffs, 0))
    return phi_arr
end

# Older 1D version with allocations
function integrate_velocities(velocities::AbstractArray{<:AbstractFloat, 1}, timediffs::Array{Int, 1})
    phi_diffs = velocities .* timediffs
    phi_arr = cumsum(coalesce.(phi_diffs, 0))
    pushfirst!(phi_arr, 0)
    return phi_arr
end

function save_line_fit(unreg_fname)
    geolist = InsarTimeseries.load_geolist_from_h5(unreg_fname)
    geolist_nums = [( g - geolist[1]).value for g in geolist]
    f = h5open(unreg_fname)
    dset = f["stack"]
    fname_out = replace(unreg_fname, ".h5" => "_linefit.h5")
    fout = h5open(fname_out, "w")
    d_create(fout, "stack", datatype(Float32), dataspace( (size(dset, 1), size(dset, 2)) ))
    dset_out = fout["stack"]
    for i in 1:size(dset, 1)
        for j in 1:size(dset, 2)
            p = polyfit(geolist_nums, vec(dset[i, j, :]), 1)
            dset_out[i, j] = p(geolist_nums[end])
        end
    end
    close(fout); close(f)
end
