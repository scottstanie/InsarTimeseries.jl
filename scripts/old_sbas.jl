"""Threaded version: If we can load all stack file into memory, might be much quicker"""
function run_sbas(unw_stack::AbstractArray{<:AbstractFloat},
                  outfile::String,
                  outdset::String,
                  geolist,
                  intlist, 
                  constant_velocity::Bool, 
                  alpha::Real,
                  L1::Bool=false,
                  prune_outliers=true,
                  sigma=3,
                  prune_fast=false;
                  demrsc=nothing,
                  reference_station=nothing,
                  ref_row=nothing,
                  ref_col=nothing,
                  window=5,
                ) 

    # Check if this file/dset already exists
    isfile(outfile) && outdset in names(outfile) && error("$outdset already exists in $outfile")

    # TODO : extract from stackavg
    if !isnothing(reference_station)
        println("Using $reference_station as reference")
        ref_row, ref_col = MapImages.station_rowcol(reference_station, demrsc)
    end
    if !(isnothing(ref_row) && isnothing(ref_col))
        println("Shifting input stack to $ref_row, $ref_col")
        # Using col, row since we 
        @time unw_stack .= InsarTimeseries.shift_stack(unw_stack, ref_row, ref_col, window=window)
    end
    L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
    nrows, ncols, _ = size(unw_stack)
 
    outstack = Array{Float32, 2}(undef, (nrows, ncols))
    println("Out size: ($nrows, $ncols)")

    @time Threads.@threads for col in 1:ncols
        for row in 1:nrows
            # soln_phase, igram_count, geo_clean, unw_clean
            unw_pixel_raw = unw_stack[row, col, :]
            soln_phase, _, _, _ = InsarTimeseries.calc_soln(unw_pixel_raw, geolist, intlist, alpha, constant_velocity;
                                            L1=L1, prune_outliers=prune_outliers, sigma=sigma,
                                            prune_fast=prune_fast)
            outstack[row, col] = InsarTimeseries.P2MM * soln_phase[1]
        end
    end

    println("Writing solution into $outfile : $outdset")
    h5write(outfile, outdset, permutedims(outstack))

    !isnothing(demrsc) && Sario.save_dem_to_h5(outfile, demrsc, overwrite=true)
    Sario.save_geolist_to_h5(outfile, outdset, geolist, overwrite=true)

    return outfile, outdset
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

# function save_line_fit(unreg_fname)
#     geolist = InsarTimeseries.load_geolist_from_h5(unreg_fname)
#     geolist_nums = [( g - geolist[1]).value for g in geolist]
#     f = h5open(unreg_fname)
#     dset = f["stack"]
#     fname_out = replace(unreg_fname, ".h5" => "_linefit.h5")
#     fout = h5open(fname_out, "w")
#     d_create(fout, "stack", datatype(Float32), dataspace( (size(dset, 1), size(dset, 2)) ))
#     dset_out = fout["stack"]
#     for i in 1:size(dset, 1)
#         for j in 1:size(dset, 2)
#             p = polyfit(geolist_nums, vec(dset[i, j, :]), 1)
#             dset_out[i, j] = p(geolist_nums[end])
#         end
#     end
#     close(fout); close(f)
# end
#
            # pix_count += 1
            # last_time = log_count(pix_count, total_pixels, 1, every=10000, last_time=last_time)
            

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


# """2nd Threaded version: If we can load all stack file into memory, might be much quicker"""
# function run_sbas(unw_stack::AbstractArray{<:AbstractFloat},
#                   outfile::String,
#                   outdset::String,
#                   geolist,
#                   intlist, 
#                   constant_velocity::Bool, 
#                   alpha::Float32,
#                   L1::Bool=false,
#                   prune::Bool=true) 
# 
#     B = prepB(geolist, intlist, constant_velocity, alpha)
#     L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
#     # TODO: do I ever really care about the abstol to change as variable?
#     rho, alpha, abstol = 1.0, 1.6, 1e-3
#     nrows, ncols, _ = size(unw_stack)
#  
#     outstack = Array{Float32, 2}(undef, (nrows, ncols))
#     countstack = similar(outstack)
#     println("Out size: ($nrows, $ncols)")
#     pix_count, total_pixels = 0, nrows * ncols
# 
#     @time Threads.@threads for col in 1:ncols
#         for row in 1:nrows
#             soln_p2mm, igram_count = calc_soln(view(unw_stack, row, col, :), B, geolist, intlist, rho, alpha, L1, prune)
#             outstack[row, col] = soln_p2mm
#             countstack[row, col] = igram_count
#         end
#     end
# 
#     println("Writing solution into $outfile : $outdset")
#     h5open(outfile, "cw") do f
#         f[outdset] = outstack
#         f[_count_dset(outdset)] = countstack
#     end
#     return outfile, outdset
# end
