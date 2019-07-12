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

