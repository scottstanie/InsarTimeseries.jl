using Dates
using Statistics
using Debugger
using LinearAlgebra

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / (-4 * π )

const STACK_FLAT_SHIFTED_DSET = "deramped_shifted_stack"


"""Runs SBAS inversion on all unwrapped igrams

Returns:

    geolist (list[datetime]): dates of each SAR acquisition

    phi_arr (ndarray): absolute phases of every pixel at each time

    deformation (ndarray): matrix of deformations at each pixel and time
"""
function run_inversion(unw_stack_file::String; outfile::String="deformation.h5", ignore_geo_file=nothing)

    if !isnothing(ignore_geo_file)
        geolist, intlist, ignore_int_indices = find_ignored(unw_stack_file, ignore_file=ignore_geo_file)
    else
        intlist = load_intlist_from_h5(unw_stack_file)
        geolist = load_geolist_from_h5(unw_stack_file)
    end

    # Prepare B matrix and timediffs used for each pixel inversion
    B = build_B_matrix(geolist, intlist)
    timediffs = day_diffs(geolist)
    if size(B, 2) != length(timediffs)
        println("Shapes of B $(size(B)) and timediffs $(size(timediffs)) not compatible")
    end


    println("Reading unw stack")
    unw_stack = load_hdf5_stack(unw_stack_file, STACK_FLAT_SHIFTED_DSET)

    vstack = invert_sbas(unw_stack, B, timediffs)
    # phi_arr = invert_sbas(unw_stack, B, timediffs)
        # geo_mask_columns=geo_mask_patch,
        # constant_vel=constant_vel,
        # alpha=alpha,
        # difference=difference,
    # )

    phi_arr = integrate_velocities(vstack, timediffs)
    # Multiply by wavelength ratio to go from phase to cm
    deformation = PHASE_TO_CM .* phi_arr

    if !isnothing(outfile)
        println("Saving deformation to $outfile")
        save_deformation(outfile, deformation, geolist)
    end
    # Now reshape all outputs that should be in stack form
    return (geolist, phi_arr, deformation)
end


"""Solve the problem Bv = d for each pixel in the stack"""
function invert_sbas(unw_stack::Array{Float32, 3}, B::Array{Float32, 2}, timediffs::Array{Int, 1})
    nrows, ncols, nlayers = size(unw_stack)
    num_geos = length(timediffs) + 1


    # Speeds up inversion to QR factorize
    # qB = qr(B)
    # Speeds up inversion to precompute pseudo inverse
    pB = pinv(B)

    # # One way: invert all columns as chunk (wont work with masking)
    # unw_cols = stack_to_cols(unw_stack)
    # velos = qB \ unw_cols
    # velos = pB * unw_cols
    # phi_cols = integrate_velocities(velos, timediffs)
    # phi_arr = cols_to_stack(phi_cols, nrows, ncols)
    

    # Pixel looping method:
    # phi_arr = Array{Float32, 3}(undef, (nrows, ncols, num_geos))

    # v = Array{Float32, 1}(undef, length(timediffs))
    vstack = Array{Float32, 3}(undef, (nrows, ncols, length(timediffs)))
    # pixel_diffs = Array{Float32, 1}(undef, size(unw_stack, 3))

    # # Column format:
    # unw_cols = copy(stack_to_cols(unw_stack))
    # vcols = Array{Float32, 2}(undef, (length(timediffs), nrows * ncols ))


    # @inbounds Threads.@threads for j in 1:size(unw_cols, 2)
    println("Using $(Threads.nthreads()) threads for invert_sbas loop")
    @inbounds Threads.@threads for j in 1:ncols
        @inbounds for i in 1:nrows
            # vstack[i, j, :] .= invert_column(unw_stack, qB, i, j)
            vstack[i, j, :] .= pB * view(unw_stack, i, j, :) 
            # vcols[:, j] .= pB * unw_cols[:, j]
            
            # phi_arr[i, j, :] .= phi_out
        end
    end
    return vstack
    # Go from velocities to phases in stack format
    # phi_arr = integrate_velocities(vstack, timediffs)

    # # Below: for column format
    # phi_cols = integrate_velocities(vcols, timediffs)
    # phi_arr = cols_to_stack(phi_cols, nrows, ncols)

    # return phi_arr
end

# function invert_column(unw_stack, pB, i::Int, j::Int)
#     c = view(unw_stack, i, j, :)
#     v = pB * c
# end

# function invert_column(pixel_diffs, qB)
#     v = qB \ pixel_diffs
# end

# open(readlines, filename)
read_geolist_file(filename::String) = sario.find_geos(filename=filename)

function find_ignored(unw_stack_file::String; ignore_file="geolist_missing.txt")
    """Read extra file to ignore certain dates of interferograms"""
    all_ints = load_intlist_from_h5(unw_stack_file)
    all_geos = load_geolist_from_h5(unw_stack_file)

    ignore_geos = sort(sario.find_geos(filename=ignore_file, parse=true))
    println("Ignoring the following .geo dates:")
    println(ignore_geos)

    ignore_igrams = [i for i in int_date_list if (i[0] in ignore_geos) || (i[1] in ignore_geos)]
    println("Ignoring $(length(ignore_igrams)) number of igrams")

    valid_geos = [g for g in geo_date_list if !(g in ignore_geos)]
    valid_igrams = [i for i in int_date_list if !(i[0] in ignore_geos) && !(i[1] in ignore_geos)]

    # ignore_geo_indices = indexin(ignore_geos, all_geos)
    ignore_int_indices = indexin(ignore_ints, all_ints)
    return valid_geos, valid_igrams, ignore_int_indices
end

function stack_to_cols(stack::Array{<:Number, 3})
    nrows, ncols, nlayers = size(stack)
    return reshape(permutedims(stack, (3, 2, 1)), (nlayers, :))
end

function cols_to_stack(columns::Array{<:Number, 2}, nrows::Int, ncols::Int)
    nlayers = size(columns, 1)
    return permutedims(reshape(columns, (nlayers, ncols, nrows)), (3, 2, 1))
end


"""Finds the number of days between successive .geo files"""
function day_diffs(geolist::Array{Date, 1})
    [difference.value for difference in diff(geolist)]
end


"""
    build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1}) -> Array{Int8, 2}

Takes the list of igram dates and builds the SBAS A matrix

Returns the incident-like matrix from the SBAS paper: A*phi = dphi
    Each row corresponds to an igram, each column to a .geo
    value will be -1 on the early (slave) igrams, +1 on later (master)
"""
function build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1})
    # We take the first .geo to be time 0, leave out of matrix
    # Match on date (not time) to find indices
    geolist = geolist[2:end]

    M = length(intlist)  # Number of igrams, 
    N = length(geolist)
    A = zeros(Int8, M, N)
    for i in 1:M
        early_date, late_date = intlist[i]
        idx_early = findfirst(isequal(early_date), geolist)
        if ~isnothing(idx_early)  # The first SLC will not be in the matrix
            A[i, idx_early] = -1
        end
        idx_late = findfirst(isequal(late_date), geolist)
        A[i, idx_late] = 1

    end
    return A
end

"""Takes the list of igram dates and builds the SBAS B (velocity coeff) matrix

Each row corresponds to an igram, each column to a .geo
Values will be t_k+1 - t_k for columns after the -1 in A,
up to and including the +1 entry
"""
function build_B_matrix(geolist::Array{Date, 1}, intlist::Array{Tuple{Date, Date}, 1})

    A = build_A_matrix(geolist, intlist)
    timediffs = day_diffs(geolist)
    return _create_B(A, timediffs)
end

function build_B_matrix(A::Array{Int8, 2}, timediffs::Array{Int, 1})
    return _create_B(A, timediffs)
end

function _create_B(A::Array{Int8, 2}, timediffs::Array{Int, 1})
    B = zeros(Float32, size(A))

    for i in 1:size(B, 1)
        row = A[i, :]
        start_idx = findfirst(row .== -1)
        # if no -1 entry, start at index 1. Otherwise
        if isnothing(start_idx)
            start_idx = 1
        else
            start_idx += 1
        end
        # +1 will always exist in row
        # End index is inclusive of the +1
        end_idx = findfirst(row .== 1)

        # Now fill in the time diffs in the index range
        B[i, start_idx:end_idx] = timediffs[start_idx:end_idx]
    end

    return B
end

"""Takes SBAS velocity output and finds phases

Arguments:
    velocities come from invert_sbas
    timediffs are the days between each SAR acquisitions
        length will be 1 less than num SAR acquisitions
"""
function integrate_velocities(velocities::AbstractArray{Float32, 1}, timediffs::Array{Int, 1})
    # multiply each column of vel array: each col is a separate solution
    phi_diffs = velocities .* timediffs

    # Now the final phase results are the cumulative sum of delta phis
    # This is equivalent to replacing missing with 0s (like np.ma.cumsum does)
    phi_arr = cumsum(coalesce.(phi_diffs, 0))

    # Add 0 as first entry of phase array to match geolist length on each col
    pushfirst!(phi_arr, 0)

    return phi_arr
end


function integrate_velocities(vstack::Array{Float32, 3}, timediffs::Array{Int, 1})
    nrows, ncols, nlayers = size(vstack)
    phi_arr = zeros(nrows, ncols, nlayers + 1)
    for j = 1:ncols
        for i = 1:nrows
            varr = view(vstack, i, j, :)
            phi_arr[i, j, :] = integrate_velocities(varr, timediffs)
        end
    end
    return phi_arr
end


"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""
function shift_stack(stack::Array{Float32, 3}, ref_row::Int, ref_col::Int; window::Int=3)
    half_win = div(window, 2)
    winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
    # temp = Array{Float32, 2}(undef, winsize)
    for k = 1:size(stack, 3)
        layer = view(stack, :, :, k)
        
        patch = layer[ref_row - half_win:ref_row + half_win, 
                      ref_col - half_win:ref_col + half_win]
        stack[:, :, k] .-= mean(patch)

    end
    return stack
end


