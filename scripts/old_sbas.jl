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
            # pix_count += 1
            # last_time = log_count(pix_count, total_pixels, 1, every=10000, last_time=last_time)
