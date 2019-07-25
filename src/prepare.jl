"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""
function shift_unw_file(unw_stack_file::String; ref_row=nothing, ref_col=nothing, window=3, ref_station=nothing, overwrite=false)
    """Runs a reference point shift on flattened stack of unw files stored in .h5"""
    if !sario.check_dset(unw_stack_file, STACK_FLAT_SHIFTED_DSET, overwrite)
        return nothing
    end

    if (isnothing(ref_row) || isnothing(ref_col))
        println("Using $ref_station as reference")
        rsc_data = sario.load_dem_from_h5(unw_stack_file)
        ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
    end
    println("Starting shift_stack: using ($ref_row, $ref_col) as reference")

    h5open(unw_stack_file, "cw") do f
        if !(STACK_FLAT_DSET in names(f))
            throw("Need $STACK_FLAT_DSET to be created in $unw_stack_file before shift stack can be run")
        end

        stack_in = f[STACK_FLAT_DSET]
        d_create(f,
            STACK_FLAT_SHIFTED_DSET,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[STACK_FLAT_SHIFTED_DSET]

        # shift_stack(stack_in, stack_out, ref_row, ref_col, window=window)
        # Note: switching these so we don't have to permute dims upon loading HDF5Dset
        shift_stack(stack_in, stack_out, ref_col, ref_row, window=window)

    end

    h5writeattr(unw_stack_file, STACK_FLAT_SHIFTED_DSET, Dict(REFERENCE_ATTR => [ref_row, ref_col]))
    h5writeattr(unw_stack_file, STACK_FLAT_SHIFTED_DSET, Dict(REFERENCE_STATION_ATTR => ref_station))

    println("Shifting stack complete")
end


"""
Note: window is size of the group around ref pixel to avg for reference.
    if window=1 or None, only the single pixel used to shift the group.
"""

# TODO: check later if worth doing one for HDF5Dataset, one for Array
function shift_stack(stack_in, stack_out, ref_row::Int, ref_col::Int; window::Int=3)
    half_win = div(window, 2)
    winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
    patch = Array{Float32, 2}(undef, winsize)

    @show stack_in
    @show stack_out
    for k = 1:size(stack_in, 3)
        println("Shiffting layer $k")
        layer = view(stack_in, :, :, k)
        
        patch .= layer[ref_row - half_win:ref_row + half_win, 
                       ref_col - half_win:ref_col + half_win]

        # Adding the `view` to eliminate extra singleton dimension from HDF5
        stack_out[:, :, k] = view(layer .- mean(patch), :, :)

    end
    return stack_out
end

function remove_ramp(z, order, mask)
    z_masked = copy(z)
    z_masked[mask] .= NaN
    return z - estimate_ramp(z_masked)

end


"""Takes a 2D array an fits a linear plane to the data
    Ignores pixels that have nan or missing values
"""
function estimate_ramp(z::Array{<:AbstractFloat, 2})
    good_idxs = .~isnan.(z)

    A = ones((sum(good_idxs), 3))
    Aidx = 1
    for idx in CartesianIndices(z)
        if good_idxs[idx]
            # row, col is equiv to y, x, but subtract 1 to start at 0
            y, x = idx.I .- 1
            A[Aidx, 2] = x
            A[Aidx, 3] = y
            Aidx += 1
        end
    end
    coeffs = A \ z[good_idxs]
    c, a, b = coeffs

    z_fit = similar(z)
    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c
    end

    return z_fit

end

function create_mean_hdf5(h5file::String; dset_name::String=STACK_DSET)
    h5open(h5file, "cw") do f
        dset = f[dset_name]
        nrows, ncols, nlayers = size(dset)
        mean_buf = zeros(eltype(dset), nrows, ncols)

        println("Creating mean from stack:")
        @show size(dset)
        for idx in 1:nlayers
            # Adding the `[:,:]` to eliminate extra singleton dimension from HDF5
            mean_buf .+= dset[:, :, idx][:, :]
        end
        mean_buf ./= nlayers
        write(f, STACK_MEAN_DSET, mean_buf)
    end
end
