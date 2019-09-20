using Distributed
using Printf
import Glob

# function run_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, 
#                   geolist,
#                   intlist, 
#                   valid_igram_indices, 
#                   constant_velocity::Bool, 
#                   alpha::Float32,
#                   L1::Bool=false) 
# 
#     B = prepB(geolist, intlist, constant_velocity, alpha)
#     @time vstack = invert_sbas(unw_stack, B, valid_igram_indices, extra_zeros=extra_zeros, L1=L1)
#     return vstack
# end

function proc_pixel(unw_stack_file, in_dset, valid_igram_indices,
                    outfile, outdset, B, geolist, intlist, rho, alpha, L1=true;
                    row=nothing, col=nothing)
    unw_pixel = h5read(unw_stack_file, in_dset, (row, col, :))[1, 1, valid_igram_indices]
    
    # Also load correlations for cutoff
    # TODO: if we dont care about corr, don't waste time loading this
    # cor_pixel = h5read(CC_FILENAME, "stack", (row, col, :))[1, 1, valid_igram_indices]
    #
    soln_p2mm, igram_count = calc_soln(unw_pixel, B, geolist, intlist, rho, alpha, L1)

    # Now with the fully clean igram list, invert
    dist_outfile = string(Distributed.myid()) * outfile
    h5open(dist_outfile, "r+") do f
        f[outdset][row, col] = soln_p2mm
        f["counts"][row, col] = igram_count
    end

end

function calc_soln(unw_pixel, B, geolist, intlist, rho, alpha, L1=true)::Tuple{Float32, Int64}
    _, intlist_clean, unw_clean, B_clean = prune_igrams(geolist, intlist, unw_pixel, B, mean_sigma_cutoff=1)

    igram_count = length(unw_clean)
    if igram_count < 50  # TODO: justify this minimum data
        soln_p2mm = Float32(0)
    else
        soln = L1 ? invert_pixel(unw_clean, B_clean, rho=rho, alpha=alpha) : B_clean \ unw_clean
        soln_p2mm = Float32.(P2MM * soln[1])
    end
    return soln_p2mm, igram_count
end

partial_files(outfile) = Glob.glob("[0-9]*" * outfile)

function merge_partial_files(outfile, dsets...)
    for dset in dsets
        tmp = zeros(size(partial_files(outfile)[1], dset))

        for f in partial_files(outfile)
            cur = h5read(f, dset)
            tmp += cur
        end
        h5open(outfile, "cw") do fout
            fout[dset] = tmp
        end
    end
    for f in partial_files(outfile)
        rm(f)
    end
end

function run_sbas(unw_stack_file::String,
                  dset::String,
                  outfile::String,
                  outdset::String,
                  geolist,
                  intlist, 
                  valid_igram_indices, 
                  constant_velocity::Bool, 
                  alpha::Float32,
                  L1::Bool=false) 

    B = prepB(geolist, intlist, constant_velocity, alpha)
    L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
    # TODO: do I ever really care about the abstol to change as variable?
    rho, alpha, abstol = 1.0, 1.6, 1e-3
    nrows, ncols, _ = size(unw_stack_file, dset)

    # h5open(outfile, "w") do fout
    #     d_create(fout, outdset, datatype(Float32), dataspace(nrows, ncols))
    # end
    println("Making new writing file for each of $(length(Distributed.workers()))")
    @time @sync @distributed for id in Distributed.workers()
        h5open(string(id)*outfile, "w") do fout
            d_create(fout, outdset, datatype(Float32), dataspace(nrows, ncols))
            d_create(fout, "counts", datatype(Int32), dataspace(nrows, ncols))
        end
    end
 
    @time @sync @distributed for (row, col) in get_unmasked_idxs()
    # @time @sync @distributed for (row, col) in collect(Iterators.product(1:1000, 1:2500))
    # @time @sync @distributed for (row, col) in collect(Iterators.product(3500:3600, 1900:2000))
        proc_pixel(unw_stack_file, dset, valid_igram_indices, outfile, 
                   outdset, B, geolist, intlist, rho, alpha, L1,
                   row=row, col=col)
    end
    println("Merging files into $outfile")
    @time merge_partial_files(outfile, outdset, "counts")
    return outfile, outdset
end

# Need the .I so we can use to load from h5 
get_unmasked_idxs(do_permute=false) = [cart_idx.I for cart_idx in findall(.!Sario.load_mask(do_permute))]

# If we can load all stack file into memory, will be much quicker
function run_sbas(unw_stack::AbstractArray{<:AbstractFloat},
                  outfile::String,
                  outdset::String,
                  geolist,
                  intlist, 
                  constant_velocity::Bool, 
                  alpha::Float32,
                  L1::Bool=false) 

    B = prepB(geolist, intlist, constant_velocity, alpha)
    L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
    # TODO: do I ever really care about the abstol to change as variable?
    rho, alpha, abstol = 1.0, 1.6, 1e-3
    nrows, ncols, _ = size(unw_stack)
 
    outstack = Array{Float32, 2}(undef, (nrows, ncols))
    countstack = similar(outstack)
    println("Out size: ($nrows, $ncols)")
    pix_count, total_pixels = 0, nrows * ncols

    # last_time = time()
    @time Threads.@threads for col in 1:ncols
        for row in 1:nrows
            soln_p2mm, igram_count = calc_soln(view(unw_stack, row, col, :), B, geolist, intlist, rho, alpha, L1)
            outstack[row, col] = soln_p2mm
            countstack[row, col] = igram_count
            # pix_count += 1
            # last_time = log_count(pix_count, total_pixels, 1, every=10000, last_time=last_time)
        end
    end

    println("Writing solution into $outfile")
    h5open(outfile, "cw") do f
        f[outdset] = outstack
        f["counts"] = countstack
    end
    return outfile, outdset
end

function prepB(geolist, intlist, constant_velocity=false, alpha=0)
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
    # TODO: solve for stacks too big for memory
    # if alpha > 0
    #     println("Regularizing solution with alpha = $alpha")
    #     B = augment_B(B, alpha)
    #     extra_zeros = size(B, 1) - length(intlist)
    # else
    #     # No regularization added to pixels
    #     extra_zeros = 0
    # end
    return B  # , extra_zeros
end

"""Solve the problem Bv = d for each pixel in the stack"""
function invert_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, B::Array{Float32, 2}, valid_igram_indices;
                     extra_zeros=0, L1=true)
    nrows, ncols, nlayers = size(unw_stack)
    total_pixels = nrows*ncols*nlayers

    num_igrams, num_timediffs = size(B)

    # Output holder
    vstack = Array{Float32, 3}(undef, (nrows, ncols, num_timediffs))

    println("Using $(Threads.nthreads()) threads for invert_sbas loop")
    # @inbounds Threads.@threads for j in 1:ncols
        # @inbounds for i in 1:nrows
    step = 1000
    row = 1
    col = 1

    # Speeds up inversion to precompute pseudo inverse for L2 least squares case
    if L1
        # v = Convex.Variable(num_timediffs)
        lu_tuple = factor(Float64.(B))
        # TODO: seems to be pretty high variance for alpha...
        # figure out which params are best/ how to adjust on fly
        # Settings found to be okay for now
        rho, alpha, abstol = 1.0, 1.6, 1e-3
    else
        pB = pinv(B)
    end

    chunk = zeros(Float32, (step, step, length(valid_igram_indices)))
    pix_count = 0
    # clear_mem_every = 10000  # TODO: figure out what causes Convex.jl to grow memory and slow down
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
            println("Reading new $row_c x $col_c chunk: ($row:$end_row, $col:$end_col)")
            @time chunk[1:row_c, 1:col_c, :] .= unw_stack[row:end_row, col:end_col, :][:, :, valid_igram_indices]

            last_time = time()
            @inbounds Threads.@threads for j in 1:col_c
                @inbounds for i in 1:row_c
                    if L1
                        # vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], B, v)
                        vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], B, rho=rho, alpha=alpha, 
                                                                    lu_tuple=lu_tuple, abstol=abstol)
                    else
                        vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], pB, extra_zeros)
                    end
                    pix_count += 1
                    last_time = log_count(pix_count, total_pixels, nlayers, every=10000, last_time=last_time)
                    # Because of a memory bug in either ECOS or Convex.jl
                    # if pix_count % clear_mem_every == 0
                        # println("CLEARING MEMORY")
                        # Convex.clearmemory()
                        # GC.gc()
                        # v = Convex.Variable(size(B, 2))
                    # end
                end
            end
            row += step
        end

        row = 1
        col += step
    end

    return vstack
end

# In case we make each pixel inversion more complicated (extra penalty functions, etc.)
# TODO: does a Val type for the case with no extra zeros help here?
function invert_pixel(pixel::Array{Float32, 1}, pB::Array{Float32,2}, extra_zeros=0)
    if extra_zeros == 0
        return pB * pixel
    else
        return pB * vcat(pixel, zeros(extra_zeros))
    end
end

# function invert_pixel(pixel::Array{Float32, 1}, B::Array{Float32,2}, v::Convex.Variable)
#     # Note: these are defined here to allow multithreading and not reuse across threads
#     solver = ECOS.ECOSSolver(verbose=0)
#     problem = Convex.minimize(norm(B*v - pixel, 1))
#     Convex.solve!(problem, solver)
# 
#     if length(v.value) > 1
#         Float32.(vec(v.value))
#     else
#         Float32.([v.value])
#     end
# end

# Defined in optimize.jl
function invert_pixel(pixel::AbstractArray{T, 1}, B::AbstractArray{T,2}; 
                      rho=1.0, alpha=1.8, lu_tuple=nothing, abstol=1e-3) where {T<:AbstractFloat}
    return Float32.(huber_fit(Float64.(B), Float64.(pixel), rho, alpha, 
                              lu_tuple=lu_tuple, abstol=abstol))
end

function log_count(pix_count, total_pixels, nlayers; every=100_000, last_time=nothing)
    if (pix_count % every) == 0
        flush(stdout)
        pct_done = 100*pix_count*nlayers/total_pixels
        @printf("Processed %.3g pixels out of %.3g (%.2f%% done).", 
                pix_count*nlayers, total_pixels, pct_done)
        if !isnothing(last_time)
            t = time()
            elapsed = t - last_time
            @printf(" %.2f seconds elapsed\n", elapsed)
            return t
        else
            @printf("\n")
        end
    else
        return last_time
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

# TODO: this can't work for huge stacks
function augment_matrices(B::Array{Float32, 2}, unw_stack::Array{Float32, 3}, alpha::Float32)
    B = Float32.(vcat(B, alpha*I))
    # Now make num rows match
    nrows, ncols, nlayers = size(unw_stack)
    zeros_shape = (nrows, ncols, size(B, 1) - nlayers)
    unw_stack = Float32.(cat(unw_stack, zeros(zeros_shape), dims=3))
    return B, unw_stack
end

