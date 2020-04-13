using SparseArrays: sparse
using LinearAlgebra: cholesky, ldiv!
# import Convex
# import ECOS


"""Fit Ax = b using Huber loss in ADMM (alternatinve directions method of multiplers)
Adapted from https://web.stanford.edu/~boyd/papers/admm/huber/huber_fit.html
"""
# TODO: differentiate this from the regularization alpha
function huber_fit(
    A,
    b,
    rho = 1.0,
    alpha = 1.0;
    lu_tuple = nothing,
    quiet = true,
    max_iter = 1000,
    abstol = 1e-4,
    reltol = 1e-2,
)
    m, n = size(A)
    Atb = A' * b

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

    kappa = 1 + 1 / rho  # For shrinkage
    r_factor = rho / (1 + rho)
    s_factor = 1 / (1 + rho)
    while idx < max_iter
        # x update
        zold .= z
        q = Atb .+ A' * (z - u)
        x .= U \ (L \ q)

        # x update with relaxation
        Ax_hat = alpha .* A * x + (1 - alpha) * (z .+ b)
        tmp = Ax_hat - b + u
        z .= (r_factor .* tmp) + (s_factor .* shrinkage(tmp, kappa))
        u += (Ax_hat - z - b)

        # Stopping check
        r_norm = norm(A * x - z - b)

        # equivalent to vv below,  but faster
        # s_norm  = norm(-rho * A' * (z - zold))
        s_norm = norm(BLAS.gemv('T', -rho, A, (z - zold)))

        eps_pri = sqrt(n) * abstol + reltol * maximum([norm(A * x), norm(-z), norm(b)])
        eps_dual = sqrt(n) * abstol + reltol * norm(rho * u)

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
    println("Caution: problem did not converge within $max_iter iterations to within $abstol, $reltol (alpha: $alpha, rho: $rho, sizes: $(size(A)), $(size(b)))")
    return x
end


"""Evaluate the huber loss with a bandwidth of `M`"""
function huber(x; M = 1)
    y = max.(x, 0)
    z = min.(y, M)
    return z .* (2 .* y .- z)
end

objective(z) = sum(huber(z)) / 2


# Sparse factoring
function factor(A)
    ch = cholesky(A' * A)
    return sparse(ch.L), sparse(ch.U)
end

pos(x) = max.(x, 0)

# Soft threshold operator (section 4.4.3)
shrinkage(x, kappa) = pos(1 .- kappa ./ abs.(x)) .* x



# ### Alternate: Solve the L1 fitting problem using Convex
# # Note: this had a memory leak earlier when run on many pixels,
# # possibly from the ECOS library
# function l1_fit(A, b, x::Union{Convex.Variable, Nothing}=nothing, quiet=true)
#     solver = ECOS.ECOSSolver(verbose=0)
#     if isnothing(x)
#         x = Convex.Variable(size(A, 2))
#     end
# 
#     problem = Convex.minimize(norm(A*x - b, 1))
#     Convex.solve!(problem, solver)
#     return x.value
# end


#### IRLS algorithm functions: #### 
l1_objective(A, x, b) = sum(abs.(A * x - b))

"""Iteratively reweighted least squares (IRLS), used to solve L1 minimization
Source: https://en.wikipedia.org/wiki/Iteratively_reweighted_least_squares"""
function irls(
    A::AbstractArray{<:AbstractFloat,2},
    b::AbstractArray{<:AbstractFloat,1};
    p::Int = 1,
    iters = 50,
)
    M, N = size(A)

    # Use Float64 to avoid roundoff NaNs
    x = Array{Float64,1}(undef, N)
    W = zeros(Float64, (size(A, 1), size(A, 1)))

    ep = sqrt(eps(eltype(x)))
    W .= diagm(0 => (abs.(b - A * x) .+ ep) .^ (p - 2))

    for ii = 1:iters
        x .= (A' * W * A) \ (A' * W * b)
        W .= diagm(0 => (abs.(b - A * x) .+ ep) .^ (p - 2))
    end
    # println("objective: ", l1_objective(A, x, b))
    return Float32.(x)
end



# # Function compatible for run_sbas
# function invert_pixel(pixel::AbstractArray{T, 1}, B::AbstractArray{T,2}; iters=50, p=1) where {T<:AbstractFloat}
#     return irls(B, pixel, iters=iters, p=p)
# end
