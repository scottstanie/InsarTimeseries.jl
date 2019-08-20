using Printf
using Convex
using ECOS
using SparseArrays: sparse
using LinearAlgebra: cholesky, ldiv!

function run_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, 
                  geolist,
                  intlist, 
                  valid_igram_indices, 
                  constant_velocity::Bool, 
                  alpha::Float32,
                  L1::Bool=false) 

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

    @time vstack = invert_sbas(unw_stack, B, valid_igram_indices, extra_zeros=extra_zeros, L1=L1)
    return vstack
end

"""Solve the problem Bv = d for each pixel in the stack"""
function invert_sbas(unw_stack::Union{HDF5Dataset, Array{Float32, 3}}, B::Array{Float32, 2}, valid_igram_indices;
                     extra_zeros=0, L1=true)
    nrows, ncols, nlayers = size(unw_stack)
    total_pixels = nrows*ncols*nlayers

    num_igrams, num_timediffs = size(B)


    # Pixel looping method:
    # TODO: maybe this should be an HDF5Dataset too?
    vstack = Array{Float32, 3}(undef, (nrows, ncols, num_timediffs))

    println("Using $(Threads.nthreads()) threads for invert_sbas loop")
    # @inbounds Threads.@threads for j in 1:ncols
        # @inbounds for i in 1:nrows
    step = 1000
    row = 1
    col = 1

    # Speeds up inversion to precompute pseudo inverse for L2 least squares case
    if L1
        # v = Variable(num_timediffs)
        lu_tuple = factor(Float64.(B))
        # Settings found to be good for these problems
        # TODO: seems to be pretty high variancce for alpha.. figure out which params are best/ how to adjust on fly
        rho, alpha, abstol = 1.0, 1.6, 1e-3
    else
        pB = pinv(B)
    end

    chunk = zeros(Float32, (step, step, length(valid_igram_indices)))
    pix_count = 0
    clear_mem_every = 10000  # TODO: figure out what causes Convex.jl to grow memory and slow down
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
                    # print("solving $i, $j")
                    if L1
                        # vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], B, v)
                        # vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], B, iters=30)
                        vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], B, rho=rho, alpha=alpha, 
                                                                    lu_tuple=lu_tuple, abstol=abstol)
                    else
                        vstack[row+i-1, col+j-1, :] .= invert_pixel(chunk[i, j, :], pB, extra_zeros)
                    end
                    pix_count += 1
                    last_time = log_count(pix_count, total_pixels, nlayers, every=10000, last_time=last_time)
                    # if pix_count % clear_mem_every == 0
                        # TODO: stupid to do this even when not using Convex...
                        # println("CLEARING MEMORY")
                        # Convex.clearmemory()
                        # GC.gc()
                        # v = Variable(size(B, 2))
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

function invert_pixel(pixel::Array{Float32, 1}, B::Array{Float32,2}, v::Convex.Variable)
    # Note: these are defined here to allow multithreading and not reuse across threads
    solver = ECOSSolver(verbose=0)
    # v = Variable(size(B, 2))

    problem = minimize(norm(B*v - pixel, 1))
    solve!(problem, solver)
    if length(v.value) > 1
        Float32.(reshape(v.value, :))
    else
        Float32.([v.value])
    end
end

function invert_pixel(pixel::AbstractArray{T, 1}, B::AbstractArray{T,2}; iters=50, p=1) where {T<:AbstractFloat}
    return irls(B, pixel, iters=iters, p=p)
end


function invert_pixel(pixel::AbstractArray{T, 1}, B::AbstractArray{T,2}; 
                      rho=1.0, alpha=1.8, lu_tuple=nothing, abstol=1e-3) where {T<:AbstractFloat}
    return Float32.(huber_fit(Float64.(B), Float64.(pixel), rho, alpha, lu_tuple=lu_tuple, abstol=abstol))
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

function augment_matrices(B::Array{Float32, 2}, unw_stack::Array{Float32, 3}, alpha::Float32)
    B = Float32.(vcat(B, alpha*I))
    # Now make num rows match
    nrows, ncols, nlayers = size(unw_stack)
    zeros_shape = (nrows, ncols, size(B, 1) - nlayers)
    unw_stack = Float32.(cat(unw_stack, zeros(zeros_shape), dims=3))
    return B, unw_stack
end

"""Iteratively reweighted least squares (IRLS), used to solve L1 minimization
Source: https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares"""
function irls(A::Array{<:AbstractFloat, 2}, b::AbstractArray{<:AbstractFloat, 1},
              W::Array{<:AbstractFloat, 2}; p::Int=1, iters=50)
    M, N = size(A)

    # Use Float64 to avoid roundoff NaNs
    x = Array{Float64, 1}(undef, N)

    _iterate_irls!(A, b, W, x, p, iters)
    # println("objective: ", l1_objective(A, x, b))
    return Float32.(x)
end
function irls(A::AbstractArray{<:AbstractFloat, 2}, b::AbstractArray{<:AbstractFloat, 1};
              p::Int=1, iters=50)
    M, N = size(A)

    # Use Float64 to avoid roundoff NaNs
    x = Array{Float64, 1}(undef, N)
    W = zeros(Float64, (size(A, 1), size(A, 1)))

    _iterate_irls!(A, b, W, x, p, iters)
    # println("objective: ", l1_objective(A, x, b))
    return Float32.(x)
end

# function _iterate_irls!(A, b, W, x, p, iters, L1, L2, Wnext)
function _iterate_irls!(A, b, W, x, p, iters)
    ep = sqrt(eps(eltype(x)))
    W .= diagm(0 => (abs.(b-A*x) .+ ep).^(p-2))

    for ii in 1:iters
        # L1 .= (A' * W * A) 
        # L2 .= (A' * W * b)
        # x .= L1 \ L2 
        # Wnext .= (abs.(b-A*x) .+ ep).^(p-2)
        # W .= diagm(0 => Wnext)
        x .= (A' * W * A) \ (A' * W * b)
        W .= diagm(0 => (abs.(b-A*x) .+ ep).^(p-2))
    end
end

l1_objective(A, x, b) = sum(abs.(A*x-b))

"""Fit Ax = b using Huber loss
Source: https://web.stanford.edu/~boyd/papers/admm/huber/huber_fit.html
"""
function huber_fit(A, b, rho=1.0, alpha=1.0; lu_tuple=nothing, quiet=true, max_iter=1000, abstol=1e-4, reltol=1e-2)
    m, n = size(A)
    Atb = A'*b

    x = zeros(n)
    z = zeros(m)
    zold = zeros(m)
    u = zeros(m)
    # If you prefactor the A matrix
    if !isnothing(lu_tuple)
        L, U = lu_tuple
    else
        L, U = factor(A)
    end

    idx = 0

    kappa = 1 + 1/rho  # For shrinkage
    r_factor = rho/(1 + rho)
    s_factor = 1/(1+rho)
    while idx < max_iter
        # x update
        zold .= z
        q = Atb .+ A' * (z - u)
        x .= U \ (L \ q)

        # x update with relaxation
        Ax_hat = alpha .* A * x + (1-alpha) * (z .+ b)
        tmp = Ax_hat - b + u
        z .= (r_factor .* tmp) + (s_factor .* shrinkage(tmp, kappa))
        u += (Ax_hat - z - b)

        # Stopping check
        r_norm = norm(A*x - z - b)

        # equivalent to vv below,  but faster
        # s_norm  = norm(-rho * A' * (z - zold))
        s_norm = norm(BLAS.gemv('T', -rho, A, (z - zold)))

        eps_pri = sqrt(n)*abstol + reltol*maximum([norm(A*x), norm(-z), norm(b)])
        eps_dual = sqrt(n)*abstol + reltol*norm(rho*u)

        if !quiet
            @show r_norm, eps_pri, s_norm, eps_dual, objective(z)
        end

        if r_norm < eps_pri && s_norm < eps_dual
            if !quiet
                println("Number of iterations to converge: $idx")
            end
            return x
        else
            idx += 1
        end

    end
    println("Caution: problem did not converge within $max_iter iterations")
    return x
end


function huber(x; M=1)
    y = max.(x, 0)
    z = min.(y, M);
    return z .* (2 .* y .- z);
end

objective(z) =  sum(huber(z)) / 2


function factor(A)
    ch = cholesky(A' * A)
    return sparse(ch.L), sparse(ch.U)
end

pos(x) = max.(x, 0)

# Soft threshold operator (section 4.4.3)
shrinkage(x, kappa) = pos(1 .- kappa ./ abs.(x)) .* x
