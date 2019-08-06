"""Subtracts reference pixel group from each layer

window is the size around the reference pixel to average at each layer
"""

function shift_unw_file(unw_stack_file::String; stack_flat_dset=nothing, order=2,
                        ref_row=nothing, ref_col=nothing, window=3, ref_station=nothing, 
                        overwrite=false)
    """Runs a reference point shift on flattened stack of unw files stored in .h5"""
    if isnothing(stack_flat_dset)
        stack_flat_dset = (order == 1) ? STACK_FLAT_DSET1 : STACK_FLAT_DSET2
    end
    stack_flat_shifted_dset = stack_flat_dset * "_shifted"

    if !sario.check_dset(unw_stack_file, stack_flat_shifted_dset, overwrite)
        return nothing
    end

    if (isnothing(ref_row) || isnothing(ref_col))
        println("Using $ref_station as reference")
        rsc_data = sario.load_dem_from_h5(unw_stack_file)
        ref_row, ref_col = gps.station_rowcol(station_name=ref_station, rsc_data=rsc_data)
    end
    println("Starting shift_stack: using ($ref_row, $ref_col) as reference")

    h5open(unw_stack_file, "cw") do f
        if !(stack_flat_dset in names(f))
            throw("Need $stack_flat_dset to be created in $unw_stack_file before shift stack can be run")
        end

        stack_in = f[stack_flat_dset]
        d_create(f,
            stack_flat_shifted_dset,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[stack_flat_shifted_dset]

        # shift_stack(stack_in, stack_out, ref_row, ref_col, window=window)
        # Note: switching these so we don't have to permute dims upon loading HDF5Dset
        shift_stack(stack_in, stack_out, ref_col, ref_row, window=window)

    end

    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_ATTR => [ref_row, ref_col]))
    # Make sure we don't write Nothing
    station_str = isnothing(ref_station) ? "" : ref_station
    h5writeattr(unw_stack_file, stack_flat_shifted_dset, Dict(REFERENCE_STATION_ATTR => station_str))

    println("Shifting stack complete")
end


"""
Note: window is size of the group around ref pixel to avg for reference.
    if window=1 or None, only the single pixel used to shift the group.
"""

# TODO: check later if worth doing one for HDF5Dataset, one for Array
function shift_stack(stack_in, stack_out, ref_row::Int, ref_col::Int; 
                     window::Int=3, chunk_layers=100)
    half_win = div(window, 2)
    winsize = (2half_win + 1, 2half_win + 1)  # Make sure it is odd sized
    patch = Array{Float32, 2}(undef, winsize)

    nrows, ncols, nlayers = size(stack_in)
    layer = Array{Float32, 2}(undef, (nrows, ncols))
    layer_out = similar(layer)

    @show stack_in
    @show stack_out

    # Read chunks in at a time to limit HDF5 reading latency

    Threads.@sync for k = 1:nlayers
        layer .= view(stack_in, :, :, k)
        Threads.@spawn layer_out .=  _shift_layer(layer, patch, ref_row, ref_col, half_win)
        stack_out[:, :, k] .= layer_out

        if k % 100 == 0
            println("Finished with $k layers")
        end
    end
    return stack_out
end

function _shift_layer(layer, patch, ref_row, ref_col, half_win)
    patch .= layer[ref_row - half_win:ref_row + half_win, 
                   ref_col - half_win:ref_col + half_win]

    # Adding the `view` to eliminate extra singleton dimension from HDF5
    return view(layer .- mean(patch), :, :)
end



"""Runs a reference point shift on flattened stack of unw files stored in .h5"""
function deramp_unw_file(unw_stack_file::String; order=2, overwrite=false, stack_flat_dset=nothing)
    if isnothing(stack_flat_dset)
        if order == 1
            stack_flat_dset = STACK_FLAT_DSET1
        elseif order == 2
            stack_flat_dset = STACK_FLAT_DSET2
        end
        println("Writing flat unws to $stack_flat_dset")
    end

    # If we already have this dataset made, skip the function
    if !sario.check_dset(unw_stack_file, stack_flat_dset, overwrite)
        return nothing
    end

    println("Starting stack deramp : using order $order")
    fmask = h5open("masks.h5", "r")
    h5open(unw_stack_file, "cw") do f
        if !(STACK_DSET in names(f))
            throw("Need $STACK_DSET to be created in $unw_stack_file before deramp stack can be run")
        end

        stack_in = f[STACK_DSET]
        nrows, ncols, nlayers = size(stack_in)

        d_create(f,
            stack_flat_dset,
            datatype(Float32),
            dataspace(size(stack_in)),
        )
        stack_out = f[stack_flat_dset]
        mask_dset = fmask[IGRAM_MASK_DSET]
    # end

        # Pre allocate buffers
        layer = similar(stack_in[:, :, 1][:, :, 1])
        layer_buf = similar(layer)
        # layer_out = similar(layer)
        mask = similar(Bool.(mask_dset[ :, :, 1][:, :, 1]))
        numel = (order == 1) ? 3 : 6
        A = ones((length(layer), numel))
        coeffs = ones(numel)
        z_fit = similar(layer)

        chunk_layers = Threads.nthreads()
        chunk =  Array{Float32, 3}(undef, (nrows, ncols, chunk_layers))
        chunk_out =  similar(chunk)
        mask_chunk =  Array{Bool, 3}(undef, (nrows, ncols, chunk_layers))
        k = 1
        while k < size(stack_in, 3)
            endk = min(k + chunk_layers - 1, lastindex(stack_in, 3))
            endc = endk - k + 1
            println("Shiffting layers $k : $endk")
            chunk[:, :, 1:endc] .= stack_in[:, :, k:endk]
            mask_chunk[:, :, 1:endc] .= Bool.(mask_dset[:, :, k:endk])

            Threads.@threads for jdx in 1:endc
                layer = chunk[ :, :, jdx]
                mask_layer = mask_chunk[:, :, jdx]
                chunk_out[:, :, jdx] .= remove_ramp(layer, order, mask_layer)

            end

            stack_out[:, :, k:endk] = chunk_out[:, :, 1:endc]
            println("sum of cunk $k is $(sum(chunk_out))")
            println("sum of stack_out $k is $(sum(stack_out[:, :, k:endk]))")

        # for k = 1:size(stack_in, 3)
            # Threads.@spawn remove_ramp(stack_out, k, layer, order, mask)
            # layer .= stack_in[ :, :, k][:, :, 1]
            # mask .= Bool.(mask_dset[:, :, k][:, :, 1])
            # @show size(layer), size(mask), size(stack_out[:,:,k])
            # Threads.@spawn remove_ramp(stack_out, k, layer, order, mask, buf=layer_buf, 
            #                                        A=A, coeffs=coeffs, z_fit=z_fit)
            #
            # wait(layer_out)
            # l = fetch(layer_out)
            # println("sum of layer $k is $(sum(l))")
            # stack_out[:, :, k] .= l

            k += chunk_layers
            if (k-1) % (1*chunk_layers) == 0
                println("Finished with $(k-1) layers")
            end
            if k > 100
                break
            end
        end

    end
    close(fmask)

    h5writeattr(unw_stack_file, stack_flat_dset, Dict("order" => order))

    println("Deramping stack complete")
end

function remove_ramp(z, order, mask; buf=nothing, A=nothing, coeffs=nothing, z_fit=nothing)
    if isnothing(buf)
        z_masked = copy(z)
    else
        z_masked = buf
        z_masked .= z
    end

    z_masked[mask] .= NaN
    if order == 1
        return z - estimate_ramp1(z_masked, A, coeffs, z_fit)
    elseif order == 2
        return z - estimate_ramp2(z_masked, A, coeffs, z_fit)
    else
        println("WARNING: Order $order not supported. Running order 1 ramp removal")
        return z - estimate_ramp1(z_masked)
    end
end


"""Takes a 2D array an fits a linear plane to the data
    Ignores pixels that have nan or missing values

For order = 1, it will be 3 numbers, a, b, c from
     ax + by + c = z
"""
# TODO: Fix the mem-allocations here
function estimate_ramp1(z::Array{<:AbstractFloat, 2}, A=nothing, coeffs=nothing, z_fit=nothing)
    good_idxs = .~isnan.(z)

    if isnothing(A)
        A = ones((length(z), 3))
    end
    if isnothing(coeffs)
        coeffs = Array{eltype(A), 1}(undef, (3,))
    end
    if isnothing(z_fit)
        z_fit = similar(z)
    end
    Aidx = 1

    for idx in CartesianIndices(z)
        if good_idxs[idx]
            # row, col is equiv to y, x, but subtract 1 to start at 0
            y, x = idx.I .- 1
            A[Aidx, 1] = x
            A[Aidx, 2] = y
            Aidx += 1
        end
    end
    coeffs .= A[good_idxs, :] \ z[good_idxs]
    a, b, c = coeffs

    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c
    end

    return z_fit
end

"""
For order = 2, it will be 6:
    a*x + b*y + c*x*y + d*x^2 + ey^2 + f
"""
# TODO: Fix the mem-allocations here
function estimate_ramp2(z::Array{<:AbstractFloat, 2}, A=nothing, coeffs=nothing, z_fit=nothing)
    good_idxs = .~isnan.(z)

    if isnothing(A)
        A = ones((length(z), 6))
    end
    if isnothing(coeffs)
        coeffs = Array{eltype(A), 1}(undef, (6,))
    end
    if isnothing(z_fit)
        z_fit = similar(z)
    end

    Aidx = 1

    for idx in CartesianIndices(z)
        if good_idxs[idx]
            # row, col is equiv to y, x, but subtract 1 to start at 0
            y, x = idx.I .- 1
            A[Aidx, 1] = x;   A[Aidx, 2] = y;  A[Aidx, 3] = x*y
            A[Aidx, 4] = x^2; A[Aidx, 5] = y^2
            Aidx += 1
        end
    end
    qA = qr(A[reshape(good_idxs, :), :])
    coeffs .=  qA \ z[good_idxs]

    a, b, c, d, e, f = coeffs

    for idx in CartesianIndices(z_fit)
        # again subtract 1 to keep consistent with the fitting
        y, x = idx.I .- 1
        z_fit[idx] = a*x + b*y + c*x*y + d*x^2 + e*y^2 + f
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
