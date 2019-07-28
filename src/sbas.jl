using Printf

function run_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, 
                  geolist,
                  intlist, 
                  valid_igram_indices, 
                  constant_velocity::Bool, 
                  alpha::Float32) 

    # Prepare A and B matrix used for each pixel inversion
    # A = build_A_matrix(geolist, intlist)
    B = build_B_matrix(geolist, intlist)
    if size(B, 2) != (length(diff(geolist)))
        println("Shapes of B $(size(B)) and geolist $(size(geolist)) not compatible")
    end

    # Only estimate 1 parameter for constant velocity B: sum all timediffs along rows
    if constant_velocity
        println("Using constant velocity for inversion solution")
        B = sum(B, dims=2)
    end
    if alpha > 0
        println("Regularizing solution with alpha = $alpha")
        B = augment_B(B, alpha)
        extra_zeros = size(B, 1) - size(unw_stack, 3)
    else
        # No regularization added to pixels
        extra_zeros = 0
    end

    @time vstack = invert_sbas(unw_stack, B, valid_igram_indices, extra_zeros=extra_zeros)
    return vstack
end

"""Solve the problem Bv = d for each pixel in the stack"""
function invert_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, B::Array{Float32, 2}, valid_igram_indices;
                     extra_zeros=0)
    nrows, ncols, nlayers = size(unw_stack)
    total_pixels = nrows*ncols*nlayers

    num_timediffs = size(B, 2)

    # Speeds up inversion to precompute pseudo inverse
    pB = pinv(B)

    # Pixel looping method:
    # TODO: maybe this should be an HDF5Dataset too?
    vstack = Array{Float32, 3}(undef, (nrows, ncols, num_timediffs))

    println("Using $(Threads.nthreads()) threads for invert_sbas loop")
    # @inbounds Threads.@threads for j in 1:ncols
        # @inbounds for i in 1:nrows
    step = 1000
    row = 1
    col = 1

    chunk = zeros(Float32, (step, step, length(valid_igram_indices)))
    pix_count = 0
    while col <= ncols
        while row <= nrows
            end_row = min(row + step - 1, nrows)
            end_col = min(col + step - 1, ncols)

            # If the chunk is not a full square, make sure we assign 
            # only part to the chunk buffer to match broadcast
            col_c = length(col:end_col)
            row_c = length(row:end_row)

            # TODO: check that reading all depth, then subselecting is actually better than
            # For-looping and reading one layer at a time
            println("Reading new chunk")
            @time chunk[1:row_c, 1:col_c, :] .= unw_stack[row:end_row, col:end_col, :][:, :, valid_igram_indices]

            for j in 1:col_c
                for i in 1:row_c
                    vstack[row+i-1, col+j-1, :] .= invert_pixel(pB, chunk[i, j, :], extra_zeros)
                    log_count!(pix_count, total_pixels, nlayers, every=100_000)
                end
            end
            row += step
        end

        row = 1
        col += step
    end

    # OLd way where we could load all in memory:
    # for j in 1:ncols
    #     for i in 1:nrows
    #         # vstack[i, j, :] .= invert_column(unw_stack, qB, i, j)
    #         vstack[i, j, :] .= pB * unw_stack[i, j, :][1, 1, :]
    #     end
    # end
    return vstack
end

# In case we make each pixel inversion more complicated (extra penalty functions, etc.)
# TODO: does a Val type for the case with no extra zeros help here?
function invert_pixel(pB, pixel::Array{Float32, 1}, extra_zeros=0)
    if extra_zeros == 0
        return pB * pixel
    else
        return pB * vcat(pixel, zeros(extra_zeros))
    end
end

function log_count!(pix_count, total_pixels, nlayers; every=100_000)
    pix_count += 1
    if (pix_count % every) == 0
        pct_done = 100*pix_count*nlayers/total_pixels
        @printf("Processed %.3g pixels out of %.3g ( %.2f %% done)\n", 
                pix_count*nlayers, total_pixels, pct_done)
    end
end

"""
    build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Igram, 1}) -> Array{Int8, 2}

Takes the list of igram dates and builds the SBAS A matrix

Returns the incident-like matrix from the SBAS paper: A*phi = dphi
    Each row corresponds to an igram, each column to a .geo
    value will be -1 on the early (slave) igrams, +1 on later (master)
"""
function build_A_matrix(geolist::Array{Date, 1}, intlist::Array{Igram, 1})
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
function build_B_matrix(geolist::Array{Date, 1}, intlist::Array{Igram, 1})

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


"""For Tikhonov regularization, pad the B matrix with alpha*I"""
function augment_B(B::Array{Float32, 2}, alpha::Float32)
    return Float32.(vcat(B, alpha*I))
end

function augment_matrices(B::Array{Float32, 2}, unw_stack::Array{Float32, 3}, alpha::Float32)
    B = Float32.(vcat(B, alpha*I))
    # Now make num rows match
    nrows, ncols, nlayers = size(unw_stack)
    zeros_shape = (nrows, ncols, size(B, 1) - nlayers)
    unw_stack = Float32.(cat(unw_stack, zeros(zeros_shape), dims=3))
    return B, unw_stack
end

