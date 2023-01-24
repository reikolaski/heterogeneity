# Compute Optimal Saving in the Income Fluctuation Problem

using .GMM, .Rouwenhorst_method, Interpolations, Plots, LaTeXStrings

default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

function EGM(params, dist=100, tol=1e-9, N=5, N_a=500, maxiter=1000)
    """
    Compute optimal policy functions using the endogenous grid method.

    Inputs
    ------
        params: parameters of the modul
        dist: initial distance to convergence
        tol: desired convergence tolerance
        N: number of income states
        N_a: number of points on the asset grid
        maxiter: maximum number of iterations

    Outputs
    -------
        c: optimal consumption policy
        a_star: optimal savings policy
        G_a: asset grid
        G_y: income grid
        Pi: income transition matrix
    """

    # Unpack parameters, discretize income space, and define some useful functions
    beta, r, gamma, a_min, rho, sigma = params
    y_grid, Pi = Rouwenhorst(rho, sigma, N)
    u(c) = c^(1 - gamma) / (1 - gamma)
    u_prime(c) = c^(-gamma)
    c_hat(B) = B.^(-1 / gamma)

    # Define (tomorrow's) asset grid and income grid. Following Kindermann and Krueger (2018), 
    # define asset grid such that there are more points near the borrowing threshold. Exponeniate 
    # residual (log) income.
    G_a  = [a_min + 2500 * ((1 + 0.08)^(i - 1) - 1) / ((1 + 0.08)^(N_a-1) - 1) for i in 1:N_a] 
    G_y = exp.(y_grid)

    # Initial consumption policy guess
    c = r * G_a .+ G_y'

    # Iterate until consumption policy function converges
    iter = 0
    a_star = zeros(length(G_a), length(G_y))
    println("\n\nRunning EGM...")
    while dist > tol && iter <= maxiter
        # RHS of the Euler equation
        B = (beta .* (1 + r) .* Pi * u_prime.(c'))'
        
        # Use budget constraint to solve for optimal savings policy
        a_star = (c_hat(B) .+ G_a .- G_y') / (1 + r)
        
        # Update consumption policy guess, either interpolating or using budget constraint (depending on asset level)
        a_star_min = a_star[1, :]
        c_old = copy(c)
        for (j, y) in enumerate(G_y)
            # Euler equation determines consumption
            interp = LinearInterpolation(a_star[:, j], c_hat(B)[:, j], extrapolation_bc=Line())
            c[:, j] = interp.(G_a)
    
            # borrowing constraint determines consumption
            c[G_a .< a_star_min[j], j] = (1 + r) .* G_a[G_a .< a_star_min[j]] .+ y .- a_min
        end

        # Update distance between consumption policy iterations
        dist = maximum(abs.(c .- c_old))

        if iter % 50 == 0
            println(String("Iteration $(iter): $(dist)"))
        end
        if iter == maxiter
            error("Did not converge")
        end
        iter += 1
    end

    return c, a_star, G_a, G_y, Pi
end

function consumption_policy_by_income()
    """
    Plot consumption as a function of wealth for someone with median labor earnings, 
    for someone with the lowest earnings, and for someone with the highest earnings.
    """

    beta = 0.96
    r = 0.02
    a_min = 0
    rho = 0.985
    sigma = sqrt(0.0346)
    gamma = 5
    params = beta, r, gamma, a_min, rho, sigma
    c, a, G_a, G_y, Pi = EGM(params)

    plot(G_a, c[:, 1], 
            label=L"y=%$(round(G_y[1], digits= 3))",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
    plot!(G_a, c[:, 3], label=L"y=%$(round(G_y[3], digits= 3))")
    plot!(G_a, c[:, 5], label=L"y=%$(round(G_y[5], digits= 3))")
    savefig("plots/consumption-policy.pdf")

    plot(G_a[G_a .< 1], c[G_a .< 1, 1], 
            label=L"y=%$(round(G_y[1], digits= 3))",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
    plot!(G_a[G_a .< 1], c[G_a .< 1, 3], label=L"y=%$(round(G_y[3], digits= 3))")
    plot!(G_a[G_a .< 1], c[G_a .< 1, 5], label=L"y=%$(round(G_y[5], digits= 3))")
    savefig("plots/consumption-policy-constrained.pdf")
end

function consumption_policy_by_CRRA_coef()
    """
    Plot consumption as a function of wealth for someone with the lowest labor earnings 
    when gamma = 2 and gamma = 10.
    """

    beta = 0.96
    r = 0.02
    a_min = 0
    rho = 0.985
    sigma = sqrt(0.0346)

    plot1 = plot()
    plot2 = plot()
    for gamma in [2, 10]
        params = beta, r, gamma, a_min, rho, sigma
        c, a, G_a, G_y, Pi = EGM(params)

        plot!(plot1, G_a, c[:, 1], 
            label=L"\gamma=%$(gamma)",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
        savefig("plots/consumption-policy-gamma.pdf")

        plot!(plot2, G_a[G_a .< 1], c[G_a .< 1, 1], 
            label=L"\gamma=%$(gamma)",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
        savefig("plots/consumption-policy-gamma-constrained.pdf")
    end
end

function consumption_policy_varying_borrowing_constraint(N=5)
    """
    Plot consumption as a function of wealth for someone with the lowest labor earnings 
    when a_min = -5 y_1, -10 y_1, and -20 y_1.
    """

    beta = 0.96
    r = 0.02
    rho = 0.985
    sigma = sqrt(0.0346)
    gamma = 5
    y_grid, Pi = Rouwenhorst(rho, sigma, N)
    y1 = exp.(y_grid)[1]

    plot1 = plot()
    plot2 = plot()
    for val in [-5, -10, -20]
        a_min = val * y1
        params = beta, r, gamma, a_min, rho, sigma
        c, a, G_a, G_y, Pi = EGM(params)

        plot!(plot1, G_a, c[:, 1],
            label=L"a_{min}=%$(val)y_1",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
        savefig("plots/consumption-policy-a_min.pdf")

        plot!(plot2, G_a[G_a .< 1], c[G_a .< 1, 1], 
            label=L"a_{min}=%$(val)y_1",
            xlabel="Assets",
            ylabel="Consumption",
            legend=:best)
        savefig("plots/consumption-policy-a_min-constrained.pdf")
    end
end

# Run experiments defined above
consumption_policy_by_income()
consumption_policy_by_CRRA_coef()
consumption_policy_varying_borrowing_constraint()