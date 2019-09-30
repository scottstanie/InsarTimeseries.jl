using Distributed
using Printf
import Glob


"""Helper to make a path to same directory as `dset`, with new name `counts`"""
_count_dset(dset) = join(["counts"; split(dset, '/')[2:end]], '/')
# _count_dset(dset) = join([split(dset, '/')[1:end-1]; "counts"], '/')
function proc_pixel(unw_stack_file, in_dset, valid_igram_indices,
                    outfile, outdset, B, geolist, intlist, rho, alpha, L1=true,
                    prune=true; row=nothing, col=nothing)
    unw_pixel = h5read(unw_stack_file, in_dset, (row, col, :))[1, 1, valid_igram_indices]
    
    # Also load correlations for cutoff
    # cor_pixel = h5read(CC_FILENAME, "stack", (row, col, :))[1, 1, valid_igram_indices]
    # cor_thresh = 0.05
    cor_pixel, cor_thresh = nothing, 0.0
    
    soln_p2mm, igram_count = calc_soln(unw_pixel, B, geolist, intlist, rho, alpha, L1;
                                       cor_pixel=cor_pixel, cor_thresh=cor_thresh, prune=prune)

    # Now with the fully clean igram list, invert
    dist_outfile = string(Distributed.myid()) * outfile
    h5open(dist_outfile, "r+") do f
        f[outdset][row, col] = soln_p2mm
        f[_count_dset(outdset)][row, col] = igram_count
    end

end

function calc_soln(unw_pixel, B, geolist, intlist, rho, alpha, L1=true;
                   cor_pixel=nothing, cor_thresh=0.0, prune=true)::Tuple{Float32, Int64}
    sigma = 3
    if !prune
        geo_clean, intlist_clean, unw_clean, B_clean = geolist, intlist, unw_pixel, B
    else
        geo_clean, intlist_clean, unw_clean, B_clean = prune_igrams(geolist, intlist, unw_pixel, B, 
                                                                    mean_sigma_cutoff=sigma,
                                                                    cor_pixel=cor_pixel,
                                                                    cor_thresh=cor_thresh)
    end

    # TODO: redo the B calculation if 
    igram_count = length(unw_clean)
    if igram_count < 50  # TODO: justify this minimum data
        soln_p2mm = Float32(0)
    else
        soln = L1 ? invert_pixel(unw_clean, B_clean, rho=rho, alpha=alpha) : B_clean \ unw_clean
        soln_p2mm = Float32.(P2MM * soln[1])
    end
    return soln_p2mm, igram_count
end

    # TODO: clean up this saving... maybe do it in a post step? handle it with config?
    if is_3d
        timediffs = day_diffs(geolist)
        println("Integrating velocities to phases")
        @time phi_arr = integrate_velocities(vstack, timediffs)
        # Multiply by wavelength ratio to go from phase to cm
        deformation = PHASE_TO_CM .* phi_arr

        println("Saving deformation to $outfile: $cur_outdset")
        @time Sario.save_hdf5_stack(outfile, cur_outdset, deformation, do_permute=true)
    end

    # TODO: I should also save the intlist... as well as the max baseline/config stuff

"""Run sbas on multiple processes"""
function run_sbas(unw_stack_file::String,
                  dset::String,
                  outfile::String,
                  outdset::String,
                  geolist,
                  intlist, 
                  valid_igram_indices, 
                  constant_velocity::Bool, 
                  alpha::Float32,
                  L1::Bool=false,
                  prune=true) 

    B = prepB(geolist, intlist, constant_velocity, alpha)
    L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
    prune ? println("Pruning .geo dates and igrams by pixel") : println("Not pruning igrams.")
    # TODO: do I ever really care about the abstol to change as variable?
    rho, alpha, abstol = 1.0, 1.6, 1e-3
    nrows, ncols, _ = size(unw_stack_file, dset)

    println("Making new writing file for each of $(length(Distributed.workers()))")
    @time @sync @distributed for id in Distributed.workers()
        h5open(string(id)*outfile, "w") do fout
            d_create(fout, outdset, datatype(Float32), dataspace(nrows, ncols))
            d_create(fout, _count_dset(outdset), datatype(Int32), dataspace(nrows, ncols))
        end
    end
 
    @time @sync @distributed for (row, col) in get_unmasked_idxs()
    # @time @sync @distributed for (row, col) in collect(Iterators.product(1:10, 1:25))
    # @time @sync @distributed for (row, col) in collect(Iterators.product(3500:3600, 1900:2000))
        proc_pixel(unw_stack_file, dset, valid_igram_indices, outfile, 
                   outdset, B, geolist, intlist, rho, alpha, L1,
                   row=row, col=col)
    end
    println("Merging files into $outfile")
    @time merge_partial_files(outfile, outdset, _count_dset(outdset))
    return outfile, outdset
end

# Need the .I so we can use to load from h5 
get_unmasked_idxs(do_permute=false) = [cart_idx.I for cart_idx in findall(.!Sario.load_mask(do_permute))]

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

# In case we make each pixel inversion more complicated (extra penalty functions, etc.)
# TODO: does a Val type for the case with no extra zeros help here?
function invert_pixel(pixel::Array{Float32, 1}, pB::Array{Float32,2}, extra_zeros=0)
    if extra_zeros == 0
        return pB * pixel
    else
        return pB * vcat(pixel, zeros(extra_zeros))
    end
end

# Defined in optimize.jl
function invert_pixel(pixel::AbstractArray{T, 1}, B::AbstractArray{T,2}; 
                      rho=1.0, alpha=1.8, lu_tuple=nothing, abstol=1e-3) where {T<:AbstractFloat}
    return Float32.(huber_fit(Float64.(B), Float64.(pixel), rho, alpha, 
                              lu_tuple=lu_tuple, abstol=abstol))
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
