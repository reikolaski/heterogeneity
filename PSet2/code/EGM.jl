using .functions, .SCF_wealth_moments, Interpolations, Plots, LaTeXStrings, DataFrames, CSV, Distributions, Random

default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

# Set parameters
beta    = 0.98
Rf      = 1.02
rho     = 0.95
sigma_y = sqrt(0.0346)
a_min   = 0
gamma   = 5
phi     = 1

# Import deterministic component of income
age_profile = DataFrame(CSV.File("data/age_profile.csv")).a
T           = length(age_profile)
age         = range(1+24, T+24)
params = phi, beta, Rf, gamma, a_min, rho, sigma_y, age_profile

# Solve problem with CRRA utility and no bequest/savings motive
c, a, G_a, G_y, Pi = EGM_finite_horizon(params)
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path, consumption_path, income_path, wealth_path, age, "plots/lifcycle-average_income-no_bequest.pdf")

# Use bisection method to calibrate phi to match average net worth at age 60 in SCF
# Note: average net worth in SCF is around 800,000---our agents' income isn't high enough to feasibly match this
model_mean_wealth_age_60 = 0
phi_low = 500000
phi_high = 1000000
phi_guess = (phi_low + phi_high) / 2
maxiter = 50
iter = 1
adjust_guess = 0.8
while abs(model_mean_wealth_age_60 - SCF_mean_wealth_age_60) > 1000 && iter <= maxiter
    # Set/update guess for phi
    global phi_guess = adjust_guess * phi_guess + (1 - adjust_guess) * ((phi_low + phi_high) / 2)

    # Solve given current guess of phi, then compute model-implied average wealth at 60
    local c, a, G_a, G_y, Pi = EGM_finite_horizon((phi_guess, beta, Rf, gamma, a_min, rho, sigma_y, age_profile))
    local wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, n_individuals=100, type="stochastic", report="mean", seed=123)
    global model_mean_wealth_age_60 = wealth_path[end]

    # If model-implied moment is lower than data moment, increase guess to save more for end-of-life; otherwise, lower guess to save less
    if model_mean_wealth_age_60 < SCF_median_wealth_age_60
        global phi_low = phi_guess
    else
        global phi_high = phi_guess
    end
    global iter += 1
end
phi = phi_guess
# Solve problem with CRRA utility and phi calibrated to match average net worth in SCF
c, a, G_a, G_y, Pi = EGM_finite_horizon((phi, beta, Rf, gamma, a_min, rho, sigma_y, age_profile))
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path[1:end-1], consumption_path[1:end-1], income_path[1:end-1], wealth_path[1:end-1], age[1:end-1], "plots/lifcycle-average-income_matching_age60_net_worth.pdf")

# Use bisection method to calibrate phi to match median net worth at age 60 in SCF
model_median_wealth_age_60 = 0
phi_low = 1
phi_high = 5000
phi_guess = (phi_low + phi_high) / 2
maxiter = 50
iter = 1
adjust_guess = 0.8
# Set/update guess for phi
while abs(model_median_wealth_age_60 - SCF_median_wealth_age_60) > 1000 && iter <= maxiter
    # Set/update guess for phi
    global phi_guess = adjust_guess * phi_guess + (1 - adjust_guess) * ((phi_low + phi_high) / 2)
    
    # Solve given current guess of phi, then compute model-implied median wealth at 60
    local c, a, G_a, G_y, Pi = EGM_finite_horizon((phi_guess, beta, Rf, gamma, a_min, rho, sigma_y, age_profile))
    local wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, n_individuals=100, type="stochastic", report="median", seed=123)
    global model_median_wealth_age_60 = wealth_path[end]

    # If model-implied moment is lower than data moment, increase guess to save more for end-of-life; otherwise, lower guess to save less
    if model_median_wealth_age_60 < SCF_median_wealth_age_60
        global phi_low = phi_guess
    else
        global phi_high = phi_guess
    end
    global iter += 1
end
phi = phi_guess
# Solve problem with CRRA utility and phi calibrated to match median net worth in SCF
c, a, G_a, G_y, Pi = EGM_finite_horizon((phi, beta, Rf, gamma, a_min, rho, sigma_y, age_profile))
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path[1:end-1], consumption_path[1:end-1], income_path[1:end-1], wealth_path[1:end-1], age[1:end-1], "plots/lifcycle-median_income-matching_age60_net_worth.pdf")

# Solve problem with Epstein-Zin utility, sigma = 1 / gamma and no bequest/savings motive
c, a, G_a, G_y, Pi = EGM_finite_horizon((1, beta, Rf, gamma, a_min, rho, sigma_y, age_profile), utility="EZ", sigma=1/gamma)
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path, consumption_path, income_path, wealth_path, age, "plots/lifcycle-average-income-EZ.pdf")

# Solve problem with Epstein-Zin utility, sigma = 2 and no bequest/savings motive
c, a, G_a, G_y, Pi = EGM_finite_horizon((1, beta, Rf, gamma, a_min, rho, sigma_y, age_profile), utility="EZ", sigma=2)
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path, consumption_path, income_path, wealth_path, age, "plots/lifcycle-average-income-EZ-phi1-sigma2.pdf")

# Solve problem with Epstein-Zin utility, sigma = 2 and phi = 5
c, a, G_a, G_y, Pi = EGM_finite_horizon((5, beta, Rf, gamma, a_min, rho, sigma_y, age_profile), utility="EZ", sigma=2)
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path[1:end-1], consumption_path[1:end-1], income_path[1:end-1], wealth_path[1:end-1], age[1:end-1], "plots/lifcycle-average-income-EZ-phi5-sigma2.pdf")

# Use bisection method to calibrate phi to match median net worth at age 60 in SCF
model_median_wealth_age_60 = 0
phi_low = 1
phi_high = 10
phi_guess = (phi_low + phi_high) / 2
maxiter = 20
iter = 1
adjust_guess = 0.8
while abs(model_median_wealth_age_60 - SCF_median_wealth_age_60) > 1000 && iter <= maxiter
    # Set/update guess for phi
    global phi_guess = adjust_guess * phi_guess + (1 - adjust_guess) * ((phi_low + phi_high) / 2)

    # Solve given current guess of phi, then compute model-implied median wealth at 60
    local c, a, G_a, G_y, Pi = EGM_finite_horizon((phi_guess, beta, Rf, gamma, a_min, rho, sigma_y, age_profile), utility="EZ", sigma=2)
    local wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, n_individuals=100, type="stochastic", report="median", seed=123)
    global model_median_wealth_age_60 = wealth_path[end]

    # If model-implied moment is lower than data moment, increase guess to save more for end-of-life; otherwise, lower guess to save less
    if model_median_wealth_age_60 < SCF_median_wealth_age_60
        global phi_low = phi_guess
    else
        global phi_high = phi_guess
    end
    global iter += 1
end
phi = phi_guess
# Solve problem with Epstein-Zin utility, sigma = 2, and phi calibrated to match median net worth in SCF
c, a, G_a, G_y, Pi = EGM_finite_horizon((phi, beta, Rf, gamma, a_min, rho, sigma_y, age_profile), utility="EZ", sigma=2)
wealth_path, consumption_path, savings_path, income_path = simulate_many_lifecycle_paths(c, a, G_a, G_y, Pi, age_profile, Rf, type="stochastic", report="mean", seed=123)
make_plot(savings_path[1:end-1], consumption_path[1:end-1], income_path[1:end-1], wealth_path[1:end-1], age[1:end-1], "plots/lifcycle-median_income-matching_age60_net_worth-EZ.pdf")
