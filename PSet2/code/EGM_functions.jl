module functions
export Rouwenhorst, EGM_finite_horizon, simulate_income_path, simulate_lifecycle_path, simulate_many_lifecycle_paths, make_plot

using Interpolations, Plots, LaTeXStrings, DataFrames, CSV, Distributions, Random

default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

function Rouwenhorst(rho, sigma, N)
    """
    Discretize an AR(1) process characterized by persistence rho and standard deviation of the
    shock sigma into N states.

    Inputs
    ------
        rho: persistence of AR(1) process
        sigma: standard deviation of shock process
        N: desired number of states

    Outputs
    -------
        y_grid: discretized grid of states
        Theta: transition matrix
    """

    # set parameters
    p = (1 + rho) / 2
    q = p
    Psi = sqrt(N - 1) * sigma

    # Initialize 2x2 matrix
    Theta_m1 = [[p 1-p];
                [1-q q]]

    # Recursively construct NxN matrix for N >= 3
    for n in 3:N
        global Theta = p .* [[Theta_m1 zeros(n-1, 1)]; [zeros(1, n-1) 0]] +
                        (1 - p) .* [[zeros(n-1, 1) Theta_m1]; [0 zeros(1, n-1)]] +
                        (1 - q) .* [[zeros(1, n-1) 0]; [Theta_m1 zeros(n-1, 1)]] +
                        q .* [[0 zeros(1, n-1)]; [zeros(n-1, 1) Theta_m1]]
        Theta[2:end-1, :] /= 2
        Theta_m1 = Theta
    end
    
    # Divide state space into grid of N discrete states
    y_grid = LinRange(-Psi, Psi, N)

    return y_grid, Theta
end

function EGM_finite_horizon(params; N=5, N_a=5000, utility::String="power", sigma=1)
    """
    Compute optimal policy functions using the endogenous grid method.

    Inputs
    ------
        params: parameters of the modul
        N: number of income states
        N_a: number of points on the asset grid
        utility: "power" or "EZ"
        sigma: elasticity of intertemporal subsitution (if utility = "EZ")

    Outputs
    -------
        c: optimal consumption policy
        a: optimal savings policy
        G_a: asset grid
        G_y: income grid
        Pi: income transition matrix
    """

    # Unpack parameters, discretize income space, and define some useful functions
    phi, beta, Rf, gamma, a_min, rho, sigma_y, age_profile = params
    y_grid, Pi = Rouwenhorst(rho, sigma_y, N)
    T = length(age_profile)

    # Define (tomorrow's) asset grid and income grid. Following Kindermann and Krueger (2018), 
    # define asset grid such that there are more points near the borrowing threshold. Note that
    # income grid is a discretization of the stochastic part of income, excluding the age 
    # dependent component.
    G_a  = [a_min + 1e6 * ((1 + 0.005)^(i - 1) - 1) / ((1 + 0.005)^(N_a-1) - 1) for i in 1:N_a] 
    G_y = exp.(y_grid)

    # Initialize policy functions
    a = zeros(length(G_a), length(G_y), T)
    c = zeros(length(G_a), length(G_y), T)

    # Consume all wealth in period T + nonworking life
    c[:, :, T] = Rf * G_a .+ age_profile[T] .* G_y'
    a[:, :, T] .= 0

    # Account for agents' need for funds at end of working life
    adj = phi

    # Power utility
    u_prime(c) = c^(-gamma)
    u_prime_inv(B) = B.^(-1 / gamma)

    # Epstein-Zin
    V = phi .* c[:, :, T]
    V_prime = phi .* Rf .* ones(length(G_a), length(G_y))

    # Solve backward given optimal behavior next period
    for t in T-1:-1:1
        c_hat = zeros(length(G_a), length(G_y))
        a_star = zeros(length(G_a), length(G_y))

        # Given a given asset level tomorrow + income shock today, compute the value of consumption that
        # satisfies the Euler equaion, then use budget constraint to solve for implied assets today.
        for (i, a) in enumerate(G_a)
            for (j, y) in enumerate(G_y)

                if utility == "power"
                    c_hat[i, j] = u_prime_inv(adj * beta * Rf * Pi[j, :]' * u_prime.(c[i, :, t+1]))
                elseif utility == "EZ"
                    c_hat[i, j] = ( beta * (Pi[j, :]' * V[i, :].^(1 - gamma))^(((sigma - 1) / sigma) / (1 - gamma) - 1) * (Pi[j, :]' * (V[i, :].^(-gamma) .* V_prime[i, :])) ).^(-sigma)
                end
                
                a_star[i, j] = (c_hat[i, j] + a - age_profile[t] * y) / Rf
            end
        end
    
        # Interpolate or use budget constraint (depending on asset level) to get consumption policy
        a_star_min = a_star[1, :]
        for (j, y) in enumerate(G_y)
            # Euler equation determines consumption
            interp = LinearInterpolation(a_star[:, j], c_hat[:, j], extrapolation_bc=Line())
            c[:, j, t] = interp.(G_a)

            # Borrowing constraint determines consumption
            c[G_a .< a_star_min[j], j, t] = Rf .* G_a[G_a .< a_star_min[j]] .+ age_profile[t] * y .- a_min

            # Optimal savings policy
            a[:, j, t] = Rf .* G_a .+ age_profile[t] * y - c[:, j, t]
        end

        # set adjustment factor to 1 for periods prior to T-1
        adj = 1

        # If using E-Z utility, update future V and V' for next period
        V_old = copy(V)
        if utility == "EZ"
            for (i, a) in enumerate(G_a)
                for (j, y) in enumerate(G_y)
                    V[i, j] = ( c[i, j, t]^((sigma - 1) / sigma) + beta * ( Pi[j, :]' * V_old[i, :].^(1 - gamma) )^(((sigma - 1) / sigma) / (1 - gamma)) ).^(sigma / (sigma - 1))
                    V_prime[i, j] = Rf * (V[i, j] / c[i, j, t])^(1 / sigma)
                end
            end
        end
    end

    return c, a, G_a, G_y, Pi
end

function simulate_income_path(G_y, Pi, age_profile; type::String="constant", init=1, seed=nothing)
    """
    Simulate lifecycle income profile.

    Inputs
    ------
        G_y: income grid
        Pi: income transition matrix
        age_profile: age dependent component of income
        type: "constant" or "stochastic"
        init: initial income states
        seed: random seed

    Outputs
    -------
        income_path: lifecycle income over T periods
        index_path: index of lifecycle income over T periods
    """

    T = length(age_profile)

    # Income path assuming constant state 
    if type == "constant"
        index_path = init * ones(Int8, T)
        income_path = G_y[index_path] .* age_profile
        return income_path, index_path

    # Stochastic income path given transition probabilities
    elseif type == "stochastic"
        index_path = zeros(Int8, T)
        i = init
        unif = rand(Random.seed!(seed), Uniform(0, 1), T)
        for t in 1:T
            index_path[t] = i
            transition_probabilities = Pi[i, :]
            cumulative_transition_probabilities = cumsum(transition_probabilities)
            
            i = findfirst(x -> x > unif[t], cumulative_transition_probabilities)
        end
        income_path = G_y[index_path] .* age_profile
        return income_path, index_path
    end
end

function simulate_lifecycle_path(c, a, G_a, G_y, Pi, age_profile, Rf; initial_assets=0, type::String="constant", init=1, seed=nothing)
    """
    Simulate lifecycle income, consumption, savings, and wealth paths.

    Inputs
    ------
        c: optimal consumption policy
        a: optimal savings policy
        G_a: asset grid
        G_y: income grid
        Pi: income transition matrix
        age_profile: age dependent component of income
        Rf: riskfree gross interest rate
        initial_assets: starting point for asset simulation
        type: "constant" or "stochastic"
        init: initial income states
        seed: random seed

    Outputs
    -------
        wealth_path: lifecycle wealth over T periods
        consumption_path: lifecycle consumption over T periods
        savings_path: lifecycle savings over T periods
        income_path: lifecycle income over T periods
    """

    # Simulate income path, then initialize vectors for wealth, consumption, and savings paths
    T = length(age_profile)
    income_path, index_path = simulate_income_path(G_y, Pi, age_profile, type=type, init=init, seed=seed)
    wealth_path = zeros(T)
    consumption_path = zeros(T)
    savings_path = zeros(T)

    # Given simulated income path, use optimal policies to determine wealth, consumption, and savings
    for t in 1:T
        wealth_path[t] = Rf * initial_assets + income_path[t]
        consumption_path[t] = c[argmin(abs.(G_a .- initial_assets)), index_path[t], t]
        savings_path[t] = a[argmin(abs.(G_a .- initial_assets)), index_path[t], t]
        initial_assets = savings_path[t]
    end

    return wealth_path, consumption_path, savings_path, income_path
end

function simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf; n_individuals=1000, initial_assets=0, type::String="stochastic", report::String="mean", seed=nothing)
    """
    Simulate lifecycle income, consumption, savings, and wealth paths for n individuals.

    Inputs
    ------
        c: optimal consumption policy
        a: optimal savings policy
        G_a: asset grid
        G_y: income grid
        Pi: income transition matrix
        age_profile: age dependent component of income
        Rf: riskfree gross interest rate
        n_individuals: number of simulations
        initial_assets: starting point for asset simulation
        type: "constant" or "stochastic"
        report: "mean" or "median" of distribution
        seed: random seed

    Outputs
    -------
        wealth_path: lifecycle wealth over T periods
        consumption_path: lifecycle consumption over T periods
        savings_path: lifecycle savings over T periods
        income_path: lifecycle income over T periods
    """

    Random.seed!(seed)

    T = length(age_profile)
    all_income_paths = zeros(n_individuals, T)
    all_wealth_paths = zeros(n_individuals, T)
    all_consumption_paths = zeros(n_individuals, T)
    all_savings_paths = zeros(n_individuals, T)

    stationary_dist = (Pi^1000)[1, :]
    cumulative_stationary_dist = cumsum(stationary_dist)
    unif_draws = rand(Uniform(0, 1), n_individuals)
    seeds = abs.(rand(Int64, n_individuals))
    for (i, unif) in enumerate(unif_draws)
        init = findfirst(x -> x > unif, cumulative_stationary_dist)
        wealth_path, consumption_path, savings_path, income_path = simulate_lifecycle_path(c, a, G_a, G_y, Pi, age_profile, Rf, initial_assets=initial_assets, type=type, init=init, seed=seeds[i])

        all_income_paths[i, :] = income_path
        all_wealth_paths[i, :] = wealth_path
        all_consumption_paths[i, :] = consumption_path
        all_savings_paths[i, :] = savings_path
    end

    if report == "mean"
        wealth_path = mean(all_wealth_paths, dims=1)'
        consumption_path = mean(all_consumption_paths, dims=1)'
        savings_path = mean(all_savings_paths, dims=1)'
        income_path = mean(all_income_paths, dims=1)'
        return wealth_path, consumption_path, savings_path, income_path
    elseif report == "median"
        wealth_path = all_wealth_paths[sortperm(all_wealth_paths[:, end]), :][floor(Int, n_individuals / 2), :]
        consumption_path = all_consumption_paths[sortperm(all_wealth_paths[:, end]), :][floor(Int, n_individuals / 2), :]
        savings_path = all_savings_paths[sortperm(all_wealth_paths[:, end]), :][floor(Int, n_individuals / 2), :]
        income_path = all_income_paths[sortperm(all_wealth_paths[:, end]), :][floor(Int, n_individuals / 2), :]
        return wealth_path, consumption_path, savings_path, income_path
    else
        return all_wealth_paths, all_consumption_paths, all_savings_paths, all_income_paths
    end
end

function make_plot(savings_path, consumption_path, income_path, wealth_path, age, savename::String)
    plot(age, savings_path / 1000, label="savings", 
        xlabel="Age",
        ylabel="Dollars (Thousands)",
        xticks=25:5:60,
        legend=:best)
    plot!(age, consumption_path / 1000, label="consumption")
    plot!(age, income_path / 1000, label="income")
    plot!(age, wealth_path / 1000, label="wealth")
    savefig(savename)
end

end