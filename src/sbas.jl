using Distributed
using StatsBase: mad
using Printf
import Glob
import Dierckx
import SparseArrays: spdiagm


"""Helper to make a group path similar to `dset`, with new base"""
_match_dset_path(dset, newgroup) = join([newgroup; split(dset, '/')[2:end]], '/')

# E.g. `counts/1` when passed `stack/1`
_count_dset(dset) = _match_dset_path(dset, "counts")
# Used if we want to track the standard deviaition 
_stddev_dset(dset) = _match_dset_path(dset, "stddev")
_stddev_raw_dset(dset) = _match_dset_path(dset, "stddev_raw")

function proc_pixel_linear(unw_stack_file, in_dset, valid_igram_indices,
                           outfile, outdset, geolist, intlist, alpha,
                           L1=true, prune=true; row=nothing, col=nothing)
    unw_pixel_raw = h5read(unw_stack_file, in_dset, (row, col, :))[1, 1, valid_igram_indices]
    if any(isnan.(unw_pixel_raw))
        println("$row, $col: nan")
        return
    end
    
    # Also load correlations for cutoff
    # cor_pixel = h5read(CC_FILENAME, "stack", (row, col, :))[1, 1, valid_igram_indices]
    # cor_thresh = 0.05
    # cor_pixel, cor_thresh = nothing, 0.0
    
    soln_phase, igram_count, geo_clean, unw_clean = calc_soln(unw_pixel_raw, geolist, intlist, alpha, true;
                                                              L1=L1, prune=prune)

    dist_outfile = string(Distributed.myid()) * outfile
    h5open(dist_outfile, "r+") do f
        f[outdset][row, col] = P2MM * soln_phase
        f[_count_dset(outdset)][row, col] = igram_count
        f[_stddev_dset(outdset)][row, col] = std(unw_clean) .* PHASE_TO_CM
        f[_stddev_raw_dset(outdset)][row, col] = std(unw_pixel_raw) .* PHASE_TO_CM
    end
    return
end

function proc_pixel_daily(unw_stack_file, in_dset, valid_igram_indices,
                          outfile, outdset, geolist, intlist, alpha,
                          L1=true, prune=true; row=nothing, col=nothing)
    unw_pixel_raw = h5read(unw_stack_file, in_dset, (row, col, :))[1, 1, valid_igram_indices]
    if any(isnan.(unw_pixel_raw))
        return
    end

    soln_velos, igram_count, geo_clean, unw_clean = calc_soln(unw_pixel_raw, geolist, intlist, alpha, false;
                                                              L1=L1, prune=prune)

    phi_arr = _unreg_to_cm(soln_velos, geo_clean, geolist)

    dist_outfile = string(Distributed.myid()) * outfile
    h5open(dist_outfile, "r+") do f
        f[outdset][row, col, :] = phi_arr
        f[_count_dset(outdset)][row, col] = igram_count
        f[_stddev_dset(outdset)][row, col] = std(unw_clean)
        f[_stddev_raw_dset(outdset)][row, col] = std(unw_pixel_raw)
    end
end

function _unreg_to_cm(soln_velos, geolist_short, geolist_full)
    # First find the phase for only the clean dates we solved for
    timediffs = day_diffs(geolist_short)
    phi_clean = integrate_velocities(soln_velos, timediffs)

    # Then linearly interpolate between these for the removed phases
    return PHASE_TO_CM * interpolate_phase(geolist_short, phi_clean, geolist_full)
end

function calc_soln(unw_pixel, geolist, intlist, alpha, constant_velocity;
                   L1=true, cor_pixel=nothing, cor_thresh=0.0, 
                   prune=true)::Tuple{Array{Float32, 1}, Int64, Array, Array}
    sigma = 3
    if !prune
        geo_clean, intlist_clean, unw_clean = geolist, intlist, unw_pixel
    else
        geo_clean, intlist_clean, unw_clean = prune_igrams(geolist, intlist, unw_pixel,
                                                           mean_sigma_cutoff=sigma,
                                                           cor_pixel=cor_pixel,
                                                           cor_thresh=cor_thresh)
    end
    igram_count = length(unw_clean)

    B = prepB(geo_clean, intlist_clean, constant_velocity, alpha)

    # Last, pad with zeros if doing Tikh. regularization
    unw_final = alpha > 0 ? augment_zeros(B, unw_clean) : unw_clean

    if igram_count < 50  # TODO: justify this minimum data
        soln_phase = [Float32(0)]
    else
        soln = L1 ? invert_pixel_l1(unw_final, B) : B \ unw_final
        soln_phase = Float32.(soln)
    end

    return soln_phase, igram_count, geo_clean, unw_clean
end


"""Run sbas on multiple processes"""
function run_sbas(unw_stack_file::String,
                  dset::String,
                  outfile::String,
                  outdset::String,
                  geolist,
                  intlist, 
                  valid_igram_indices, 
                  constant_velocity::Bool, 
                  alpha::Real,
                  L1::Bool=false,
                  prune=true) 

    L1 ? println("Using L1 penalty for fitting") : println("Using least squares for fitting")
    prune ? println("Pruning .geo dates and igrams by pixel") : println("Not pruning igrams.")
    alpha > 0 ? println("Regularizing solution with alpha = $alpha") : println("No regularization")

    # TODO: do I ever really care about these to change as variable?
    # rho, alpha, abstol = 1.0, 1.6, 1e-3
    nrows, ncols, _ = size(unw_stack_file, dset)

    # Choose the correct processing function and size based on if we choose constant
    # TODO: is using multiple dispatch here better in any way (clarity/speed)?
    if constant_velocity
        println("Using constant velocity for inversion solution")
        proc_func = proc_pixel_linear
        outsize = (nrows, ncols)
    else
        println("Not using linear: Finding solution for each day")
        proc_func = proc_pixel_daily
        outsize = (nrows, ncols, length(geolist))
    end

    println("Making new writing file for each of $(length(Distributed.workers()))")
    @time @sync @distributed for id in Distributed.workers()
        h5open(string(id)*outfile, "w") do fout
            d_create(fout, outdset, datatype(Float32), dataspace(outsize))
            d_create(fout, _count_dset(outdset), datatype(Int32), dataspace(nrows, ncols))
            d_create(fout, _stddev_dset(outdset), datatype(Float64), dataspace(nrows, ncols))
            d_create(fout, _stddev_raw_dset(outdset), datatype(Float64), dataspace(nrows, ncols))
        end
    end
 
    @time @sync @distributed for (row, col) in get_unmasked_idxs(geolist)
    # @time @sync @distributed for (row, col) in collect(Iterators.product(195:200, 195:200))
        proc_func(unw_stack_file, dset, valid_igram_indices, outfile, 
                   outdset, geolist, intlist, alpha, L1, prune,
                   row=row, col=col)
    end
    println("Merging files into $outfile")
    @time merge_partial_files(outfile, outdset, _count_dset(outdset), 
                              _stddev_dset(outdset), _stddev_raw_dset(outdset))

    _save_std_attrs(unw_stack_file, outdset)
    return outfile, outdset
end

# Need the .I so we can use to load from h5 
get_unmasked_idxs(geolist, do_permute=false) = [cart_idx.I for cart_idx in
                                                findall(.!Sario.load_mask(geolist, do_permute=do_permute))]

function _save_std_attrs(unw_stack_file, outdset)
    attr_dict1 = Dict("description" => "original unwrapped (shifted) std. dev of pixels", "units" => "cm")
    h5writeattr(unw_stack_file, _stddev_raw_dset(outdset), attr_dict1)
    attr_dict2 = Dict("description" => "Final, outlier-removed std. dev of pixels", "units" => "cm")
    h5writeattr(unw_stack_file, _stddev_dset(outdset), attr_dict2)
end

"""Linearly interpolate the values between dates that we didn't solve for"""
function interpolate_phase(geo_clean, phis, geolist_full)
    return interp1d(_get_day_nums(geo_clean), phis, _get_day_nums(geolist_full))
end

function interp1d(x, y, xv, k=1)
    sp = Dierckx.Spline1D(x, y, k=k)
    return sp(xv)
end


function prepB(geolist, intlist, constant_velocity=false, alpha=0)
    if constant_velocity
        return Float32.(reshape(temporal_baseline(intlist), :, 1))
    end
    # Prepare A and B matrix used for each pixel inversion
    # A = build_A_matrix(geolist, intlist)
    B = build_B_matrix(geolist, intlist)

    if alpha > 0
        B = augment_B(B, alpha)
    end
    return B
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
function invert_pixel_l1(pixel::AbstractArray{<:AbstractFloat}, B::AbstractArray{<:AbstractFloat}; 
                      rho=1.0, alpha=1.6, lu_tuple=nothing, abstol=1e-3)
    return Float32.(huber_fit(Float64.(B), Float64.(pixel), rho, alpha, 
                              lu_tuple=lu_tuple, abstol=abstol))
end


"""Takes the list of igram dates and builds the SBAS B (velocity coeff) matrix

Each row corresponds to an igram, each column to a .geo
Values will be t_k+1 - t_k for columns after the -1 in A,
up to and including the +1 entry
"""
function build_B_matrix(geolist::Array{Date, 1}, intlist::Array{Igram, 1})

    timediffs = day_diffs(geolist)
    # We take the first .geo to be time 0, leave out of matrix
    # Match on date (not time) to find indices
    geolist2 = @view geolist[2:end]

    M = length(intlist)  # Number of igrams, 
    N = length(geolist2)

    B = zeros(Float32, (M, N))

    for i = 1:M
        # row = A[i, :]
        early_date, late_date = intlist[i]
        idx_early = findfirst(isequal(early_date), geolist2)
        if isnothing(idx_early)  # The first SLC will not be in the matrix
            start_idx = 1
        else
            start_idx = idx_early + 1
        end
        idx_late = findfirst(isequal(late_date), geolist2)
        end_idx = idx_late

        # Now fill in the time diffs in the index range
        B[i, start_idx:end_idx] .= @view timediffs[start_idx:end_idx]
    end

    return B
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



_D(n) = spdiagm(0 => -ones(n-1), 1=>ones(n-1))[1:end-1, :]

"""For Tikhonov regularization, pad the B matrix with alpha*D"""
function augment_B(B::Array{Float32, 2}, alpha::Real)
    extra = alpha*I
    n = size(B, 2)
    return Float64.(vcat(B, alpha * _D(n)))  # Note: as of 9/29/19, qr has bug for sparse float32 matrices
end

augment_zeros(B, unw) = vcat(unw, zeros(size(B, 1) - length(unw)))

# function augment_B(B::Array{Float32, 2}, alpha::Float32)
#     extra = alpha*I
#     return Float32.(vcat(B, ))
# end
 
# TODO: this can't work for huge stacks
function augment_matrices(B::Array{Float32, 2}, unw_stack::Array{Float32, 3}, alpha::Real)
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
function integrate_velocities(velocities::AbstractArray{<:AbstractFloat, 1}, timediffs::Array{Int, 1})
    phi_diffs = velocities .* timediffs
    phi_arr = cumsum(coalesce.(phi_diffs, 0))
    pushfirst!(phi_arr, 0)
    return phi_arr
end


# TODO: I should also save the intlist... as well as the max baseline/config stuff

### Outlier functions ######
# TODO: figure out which should go in separate file for cleanliness
function prune_igrams(geolist, intlist, unw_pixel;  # B;
                      cor_pixel=nothing, cor_thresh=0.0,
                      mean_sigma_cutoff=3, 
                      fast_cm_cutoff=SENTINEL_WAVELENGTH/4)
    # 3 Reasons for removing igrams:
    # 1. remove all from .geo dates with huge averages (whole day is garbage)
    # 2. remove very low correlation
    #   NOTE: for now, we are skipping... doesn't seem to make much difference generally
    # 3. remove longer igrams w/ longer baseline than is possible
    #   e.g. if we have 2.5 cm/year velocity, anything longer than ~1 year
    #   will be all noise- we can't reliably sense such quickly moving ground
    # TODO: verify this third criteria?

    # @show std(unw_pixel)
    # @show "start", size(intlist)
    # 1. find outliers in this pixels' values and remove them
    geo_clean, intlist_clean, unw_clean = peel_nsigma(geolist, intlist, unw_pixel, nsigma=mean_sigma_cutoff)
    # geo_clean, intlist_clean, unw_clean, B_clean = peel_nsigma(geolist, intlist, unw_pixel, B, nsigma=mean_sigma_cutoff)

    # @show "outlier: ", size(intlist_clean)

    # 2. low correlation cleaning
    if cor_thresh > 0 && !isnothing(cor_pixel)
        low_cor_igrams = intlist[cor_pixel .< cor_thresh]
        intlist_clean, unw_clean = remove_igrams(intlist_clean, unw_clean, low_cor_igrams)
    end
    # @show "cor: ", size(intlist_clean)

    # 3. with rought velocity estimate, find igrams with too long of baseline
    # Here we assume beyond some wavelength fraction is too long to sense in one igram
    if fast_cm_cutoff < 5
        Blin = prepB(geo_clean, intlist_clean, true)
        velo_cm = PHASE_TO_CM * (Blin \ unw_clean)[1]  # cm / day
        # rho, alpha, abstol = 1.0, 1.6, 1e-3
        # velo_cm = PHASE_TO_CM * invert_pixel_l1(unw_clean, B_clean, rho=rho, alpha=alpha, abstol=abstol)[1]
        day_cutoff = fast_cm_cutoff / abs(velo_cm)
        # @show (365*velo_cm), fast_cm_cutoff, day_cutoff
        too_long_igrams = [ig for ig in intlist_clean 
                           if temporal_baseline(ig) > day_cutoff]

        intlist_clean, unw_clean = remove_igrams(intlist_clean, unw_clean, too_long_igrams)
        # intlist_clean, unw_clean, B_clean  = remove_igrams(intlist_clean, unw_clean, B_clean, too_long_igrams)
    end
    # @show "fast: ", size(intlist_clean)
    # @show std(unw_clean)
    return geo_clean, intlist_clean, unw_clean
end

function remove_igrams(intlist, unw_vals,
                       bad_items::Union{Date, AbstractArray{Date}, AbstractArray{Igram}}...)
    good_idxs = trues(size(intlist))
    for item in bad_items
        good_idxs .&= _good_idxs(item, intlist)
    end
    return remove_by_idx(good_idxs, intlist, unw_vals)
end

"""Pass directly the indexes of good ones igrams (removes all BUT the good_idxs)"""
function remove_by_idx(good_idxs::AbstractArray, intlist, unw_vals)
    unw_clean = unw_vals[good_idxs]
    intlist_clean = intlist[good_idxs]
    return intlist_clean, unw_clean
end

# If we pass the B, also remove rows for bad idxs
function remove_by_idx(good_idxs::AbstractArray, intlist, unw_vals, B)
    B_clean = B[good_idxs, :]
    unw_clean = unw_vals[good_idxs]
    intlist_clean = intlist[good_idxs]
    return intlist_clean, unw_clean, B_clean 
end

function _good_idxs(bad_date_arr::AbstractArray{Date}, intlist)
    return .!reduce((x, y)-> x .| y, 
                    (b in intlist for b in bad_date_arr),
                    init=falses(size(intlist)));
end

function _good_idxs(bad_igram_arr::AbstractArray{Igram}, intlist)
    return [!(ig in bad_igram_arr) for ig in intlist]
end

_good_idxs(bad_date::Date, intlist) = .!(bad_date in intlist)



""" 3* the sigma valud as calculated using MAD
Used for robust variance est. (as a cutoff for outliers)
See https://en.wikipedia.org/wiki/Robust_measures_of_scale#IQR_and_MAD"""
mednsigma(arr, n=3) = n * mad(arr, normalize=true)  

two_way_cutoff(arr, nsigma) = (min(0, median(arr) - mednsigma(arr, nsigma)),  # Make sure it's negative
                               median(arr) + mednsigma(arr, nsigma))


"""Simpler means by day: just abs. value of phases, either mean or median"""
mean_abs_val(geolist, intlist, unw_vals) = [mean(abs.(vals_by_date(d, intlist, unw_vals)))
                                            for d in geolist];
median_abs_val(geolist, intlist, unw_vals) = [median(abs.(vals_by_date(d, intlist, unw_vals)))
                                              for d in geolist];

"""Get the values (unw, cc, etc.) of one date"""
vals_by_date(date::Date, intlist, vals) = vals[date in intlist]

"""Get the values (unw, cc, etc.) grouped by each date"""
vals_by_date(date_arr::Array{Date}, intlist, vals) = [vals[d in intlist] for d in date_arr]

function two_way_outliers(arr, nsigma)
    low, high = two_way_cutoff(arr, nsigma)
    return (arr .< low) .| (arr .> high)
end


"""Remove all igrams corresponding to the dates with means falling outside an `nsigma` interval"""
function peel_nsigma(geo, int, val; nsigma=3)
    dates_to_remove = nsigma_days(geo, int, val, nsigma)
    int2, val2 = remove_igrams(int, val, dates_to_remove)
    geo2 = [g for g in geo if !(g in dates_to_remove)]
    return geo2, int2, val2
end

"""Return the days of `geo` which are more than nsigma away from mean"""
function nsigma_days(geo, int, val, nsigma=3)
    means = mean_abs_val(geo, int, val)
    # means = median_abs_val(geo, int, val)  #; println("median abs")
    # @show mednsigma(means, nsigma)
    # TODO: change "means" if its not means
    # means = abs.(oneway_val(geo, int, val, mean))
    # means = abs.(oneway_val(geo, int, val, median))
    # means = oneway_val(geo, int, val, mean)
    
    # means = oneway_val(geo, int, val, median); # println("median oneway")
    # @show mednsigma(means, nsigma)
    # TODO: change "means" if its not means
    out_idxs = two_way_outliers(means, nsigma)

    # # FOR PRINTING ONLY
    # low, high = two_way_cutoff(means, nsigma)
    # println("Using cutoff around $(median(means)), spread $(mednsigma(means, nsigma)) : ($low, $high) ")
    # println("Number removed: $(sum(out_idxs))")
    return geo[out_idxs]
end
