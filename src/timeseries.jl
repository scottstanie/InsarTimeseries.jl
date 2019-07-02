using Dates
using Statistics
using Debugger
using LinearAlgebra

const SENTINEL_WAVELENGTH = 5.5465763  # cm
const PHASE_TO_CM = SENTINEL_WAVELENGTH / (-4 * π )


"""Runs SBAS inversion on all unwrapped igrams

Returns:
    geolist (list[datetime]): dates of each SAR acquisition from read_geolist
    phi_arr (ndarray): absolute phases of every pixel at each time
    deformation (ndarray): matrix of deformations at each pixel and time
"""
function run_inversion(igram_path::String, reference::Tuple{Int, Int}; window::Int=3, unw_ext::String=".unwflat")

    intlist = sario.read_intlist(filepath=igram_path)
    geolist = sario.read_geolist(filepath=igram_path)

    # Prepare B matrix and timediffs used for each pixel inversion
    B = build_B_matrix(geolist, intlist)
    timediffs = day_diffs(geolist)
    if size(B, 2) != length(timediffs)
        println("Shapes of B $(size(B)) and timediffs $(size(timediffs)) not compatible")
    end

    println("Reading unw stack")
    unw_stack = load_stack(directory=igram_path, file_ext=unw_ext)

    ref_row, ref_col = reference

    println("Starting shift_stack: using $ref_row, $ref_col as ref_row, ref_col")
    @time unw_stack = shift_stack(unw_stack, ref_row, ref_col, window=window)
    println("Shifting stack complete")

    phi_arr = invert_sbas(unw_stack, B, timediffs)
        # geo_mask_columns=geo_mask_patch,
        # constant_vel=constant_vel,
        # alpha=alpha,
        # difference=difference,
    # )

    # Multiple by wavelength ratio to go from phase to cm
    deformation = similar(phi_arr)
    @. deformation = PHASE_TO_CM * phi_arr

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
    phi_arr = Array{Float32, 3}(undef, (nrows, ncols, num_geos))
    # v = Array{Float32, 1}(undef, length(timediffs))
    vstack = Array{Float32, 3}(undef, (nrows, ncols, length(timediffs)))
    # pixel_diffs = Array{Float32, 1}(undef, size(unw_stack, 3))

    # # Column format:
    # unw_cols = copy(stack_to_cols(unw_stack))
    # vcols = Array{Float32, 2}(undef, (length(timediffs), nrows * ncols ))


    # @inbounds Threads.@threads for j in 1:size(unw_cols, 2)
    println("Using $(Threads.nthreads()) nubmer of threads")
    @inbounds Threads.@threads for j in 1:ncols
        @inbounds for i in 1:nrows
            # vstack[i, j, :] .= invert_column(unw_stack, qB, i, j)
            vstack[i, j, :] .= pB * view(unw_stack, i, j, :) 
            # vcols[:, j] .= pB * unw_cols[:, j]
            
            # phi_arr[i, j, :] .= phi_out
        end
    end
    # Go from velocities to phases in stack format
    phi_arr = integrate_velocities(vstack, timediffs)
    # # Below: for column format
    # phi_cols = integrate_velocities(vcols, timediffs)
    # phi_arr = cols_to_stack(phi_cols, nrows, ncols)

    return phi_arr
end

# function invert_column(unw_stack, pB, i::Int, j::Int)
#     c = view(unw_stack, i, j, :)
#     v = pB * c
# end

# function invert_column(pixel_diffs, qB)
#     v = qB \ pixel_diffs
# end


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
function integrate_velocities(velocities::Array{Float32, 1}, timediffs::Array{Int, 1})
    # multiply each column of vel array: each col is a separate solution
    phi_diffs = similar(velocities)
    phi_diffs .= velocities .* timediffs

    # Now the final phase results are the cumulative sum of delta phis
    # This is equivalent to replacing missing with 0s (like np.ma.cumsum does)
    phi_arr = cumsum(coalesce.(phi_diffs, 0))

    # Add 0 as first entry of phase array to match geolist length on each col
    pushfirst!(phi_arr, 0)

    return phi_arr
end

function integrate_velocities!(phi_arr::Array{Float32, 1}, phi_diffs::Array{Float32, 1}, phi_out::Array{Float32, 1}, velocities::Array{Float32, 1}, timediffs::Array{Int, 1})
    # multiply each column of vel array: each col is a separate solution
    # phi_diffs = similar(velocities)
    phi_diffs .= velocities .* timediffs

    # Now the final phase results are the cumulative sum of delta phis
    # This is equivalent to replacing missing with 0s (like np.ma.cumsum does)
    cumsum!(phi_arr, coalesce.(phi_diffs, 0))

    # Add 0 as first entry of phase array to match geolist length on each col
    # phi_out[1] = 0
    for i in 1:length(phi_arr)
        phi_out[i+1] = phi_arr[i]
    end
    # pushfirst!(phi_arr, 0)

    return phi_out
end

function integrate_velocities(vstack::Array{Float32, 3}, timediffs::Array{Int, 1})
    # multiply each depth stack of vel array
    phi_diffs = vstack .* timediffs[1, 1, :]

    phi_arr = cumsum(coalesce.(phi_diffs, 0), dims=3)

    # Add 0 as first entry of phase array to match geolist length on each col
    return cat(zeros(eltype(phi_arr), (size(vstack, 1), size(vstack, 2))), phi_arr, dims=3)
end

# Separate one for 2D velocities array
function integrate_velocities(velocities::Array{Float32, 2}, timediffs::Array{Int, 1})
    nrows, ncols = size(velocities)
    phi_arr = Array{Float32, 2}(undef, (nrows + 1, ncols))    
    col = Array{Float32, 1}(undef, (nrows,))
    for j in 1:ncols
        col .= velocities[:, j]
        phi_arr[:, j] .= integrate_velocities(col, timediffs)
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


