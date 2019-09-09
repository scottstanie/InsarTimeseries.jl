using Dates
using Statistics
using LinearAlgebra

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / (-4 * π )
const P2MM = 365 * 10 * PHASE_TO_CM   # mm / year


"""Runs SBAS inversion on all unwrapped igrams

Returns:

    geolist (list[datetime]): dates of each SAR acquisition

    phi_arr (ndarray): absolute phases of every pixel at each time

    deformation (ndarray): matrix of deformations at each pixel and time
"""
unw_stack
function run_inversion(unw_stack_file::String; 
                       outfile::Union{String,Nothing}=nothing, 
                       use_stackavg::Bool=false, 
                       constant_velocity::Bool=true, 
                       ignore_geo_file=nothing, 
                       max_temporal_baseline::Union{Int, Nothing}=nothing,
                       alpha::Float32=0.0f0,
                       L1::Bool=false)

    # the valid igram indices is out of all layers in the stack and mask files 
    geolist, intlist, valid_igram_indices = load_geolist_intlist(unw_stack_file, ignore_geo_file, max_temporal_baseline)


    if use_stackavg
        is_hdf5 = false
        flat_dset = STACK_FLAT_DSET
        println("Averaging stack for solution")
        vstack = run_stackavg(unw_stack_file, flat_dset, geolist, intlist)
        outdset = "stack"
        is_3d = true  # TODO: stack avg should really be 2d velo
    else
        is_hdf5 = true
        flat_dset = STACK_FLAT_SHIFTED_DSET
        println("Performing SBAS solution")
        println("Reading unw stack")
        # @time unw_stack = load_hdf5_stack(unw_stack_file, STACK_FLAT_SHIFTED_DSET)
        # @time unw_stack = load_hdf5_stack(unw_stack_file, STACK_FLAT_SHIFTED_DSET, valid_igram_indices)
        # TODO: also do this for the masks

        # vstack = run_sbas(unw_stack, geolist, intlist, constant_velocity, alpha)
        # outdset = "stack"
        # h5open(unw_stack_file) do unw_file
        #     unw_stack = unw_file[flat_dset]
        #     vstack = run_sbas(unw_stack, geolist, intlist, valid_igram_indices,
        #                       constant_velocity, alpha, L1)
        # end
        outdset = "velos"
        velo_file_out = run_sbas(unw_stack_file, flat_dset, outfile, outdset,
                                 geolist, intlist, valid_igram_indices, 
                                 constant_velocity, alpha, L1)
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
    if !isnothing(outfile)
        if is_3d
            println("Saving deformation to $outfile")
            @time save_deformation(outfile, deformation, dem_rsc, unw_stack_file=unw_stack_file, do_permute=!is_hdf5)
        else
            sario.save_dem_to_h5(outfile, dem_rsc, dset_name=DEM_RSC_DSET, overwrite=true)
        end

        save_geolist_to_h5(outfile, geolist, overwrite=true)
        save_reference(outfile, unw_stack_file, outdset, flat_dset)
    end

    # return (geolist, phi_arr, deformation)
end


function load_geolist_intlist(unw_stack_file, ignore_geo_file, max_temporal_baseline)
    geolist = load_geolist_from_h5(unw_stack_file)
    intlist = load_intlist_from_h5(unw_stack_file)

    # If we are ignoreing some indices, remove them from for all pixels
    valid_geo_indices, valid_igram_indices = find_valid_indices(geolist, intlist, ignore_geo_file, max_temporal_baseline)
    return geolist[valid_geo_indices], intlist[valid_igram_indices], valid_igram_indices
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
