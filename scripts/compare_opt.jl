import SCS
import Convex
import InsarTimeseries

import PyPlot
plt = PyPlot

huber_fit = InsarTimeseries.huber_fit
huber = InsarTimeseries.huber

obj_l1(z) = norm(z, 1)
obj_l2(z) = norm(z, 2)
obj_huber(z) = sum(huber(z)) / 2

function l1_fit(A, b)
    x = Convex.Variable(size(A, 2))
    problem = Convex.minimize(norm(A*x - b, 1))
    Convex.solve!(problem, SCS.SCSSolver(verbose=false))
    return x.value
end

# Fit a line with L1, huber, and L2
function compare_l1_huber(m=500, m_true=1, b_true=-5, sigma=3; with_outlier=true)
    # slope of line is m_true, y intercept is b_true
    x_true = collect(range(0, 10, length=m))
    data_true = (m_true .* x_true) .+ b_true
    @show m_true, b_true
    # b = x_true .+ (sigma .* randn(m))  # Gaussian noise
    b = data_true .+ (2sigma .* rand(m) .- sigma)  # Uniform noise [-sigma, sigma]
    if with_outlier
        extra = maximum(b) - minimum(b)
        b[1] += 2extra
        b[end] -= 2extra
    end

    n = 1
    A =  hcat(x_true, ones(m))

    huber_sol = huber_fit(A, b)
    l1_sol = l1_fit(A, b)
    l2_sol = A \ b

    plt.figure()
    plt.plot(x_true, data_true, label="truth")
    plt.scatter(x_true, b, label="data")
    plt.plot(x_true, A*huber_sol, label="Huber")
    plt.plot(x_true, A*l1_sol, label="L1")
    plt.plot(x_true, A*l2_sol, label="L2")
    plt.legend()
    plt.show(block=false)

    res_huber = A*huber_sol .- b
    res_l1    = A*l1_sol .- b
    res_l2    = A*l2_sol .- b
    return huber_sol, l1_sol, l2_sol, huber_sol - l1_sol, res_huber, res_l1, res_l2
end


function plot_resid(rh, r1, r2, bins=nothing)
    vm = max(maximum(abs.(rh)), maximum(abs.(r1)))
    fig, axes = plt.subplots(1, 3)
    n1, _, _ = axes[1].hist(rh, range=(-vm, vm), bins=bins)
    n2, _, _ = axes[2].hist(r1, range=(-vm, vm), bins=bins)
    n3, _, _ = axes[3].hist(r2,  range=(-vm, vm), bins=bins)
    axes[1].set_title("Huber residuals")
    axes[2].set_title("L1 (CVX) residuals")
    axes[3].set_title("L2 least sq.")
    nmax = max(maximum(n1), maximum(n2), maximum(n3))
    axes[1].set_ylim(0, nmax+2)
    axes[2].set_ylim(0, nmax+2)
    axes[3].set_ylim(0, nmax+2)
    return fig, axes
end

points = 40
sh, s1, s2, dif, rh, r1, r2 = compare_l1_huber(points, with_outlier=true)

nbins = 35
fig, axes = plot_resid(rh, r1, r2, nbins)
