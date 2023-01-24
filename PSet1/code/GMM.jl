# Estimate idiosyncratic labor income process parameters using GMM,
# to minimize distance between empirical and theoretical variances and
# autocovariances.

module GMM
export rho_GMM, var_alpha_GMM, var_nu_GMM, var_eta_GMM
using CSV, DataFrames, StatFiles, Statistics, LinearAlgebra, Optim

# Use cleaned data file exported by clean_data.do
reslogwage = DataFrame(CSV.File("data/residual_logwage.csv"))

function empirical_moments(data)
    """
    Compute empirical moments from data.

    Inputs
    ------
        data: observations are individuals with each column their wage at ages j=0,...,J

    Outputs
    -------
        moments: vector of nonredundant variance and covariance moments
    """

    T = ncol(data)
    cov_mat = zeros(T, T)

    # compute variances and all autocovariances using observations where both variables are non-missing
    for (i, c1) in enumerate(eachcol(data))
        for (j, c2) in enumerate(eachcol(data))
            sx, sy = skipmissings(c1, c2)
            cov_mat[i, j] = cov(collect(sx), collect(sy))
        end
    end

    # keep upper triangular of the covariance matrix and stack into a vector
    moments = cov_mat[triu!(trues(size(cov_mat)))], T
    return moments
end

function theoretical_moments(params, data)
    """
    Compute theoretical moments using given parameters.

    Inputs
    ------
        params: parameters of the stochastic income process
        data: observations of J wages per individual

    Outputs
    -------
        moments: vector of nonredundant variance and covariance moments
    """

    rho, var_alpha, var_nu, var_eta = params
    T = ncol(data)
    cov_mat = zeros(T, T)

    # Looping over pairs of dates, compute moments according to formulae (see writeup for details).
    for t in 0:T-1
        if t == 0
            cov_mat[t+1, t+1] = var_alpha + var_nu
        else
            cov_mat[t+1, t+1] = var_alpha + var_nu + var_eta * sum([rho^(2*(t - k)) for k in 1:t])
        end
        for th in t+1:T-1
            if t == 0
                cov_mat[t+1, th+1] = var_alpha
            else
                cov_mat[t+1, th+1] = var_alpha + var_eta * rho^(th - t) * sum([rho^(2*(t - k)) for k in 1:t])
            end
        end
    end

    # keep upper triangular of the covariance matrix and stack into a vector
    moments = cov_mat[triu!(trues(size(cov_mat)))]
    return moments
end

function GMM_objective(params, Mhat, data, W)
    """
    Compute weighted distance between theoretical and empirical moments.

    Inputs
    ------
        params: parameters of the stochastic income process
        Mhat: empirical moments
        data: observations of J wages per individual
        W: weighting matrix

    Outputs
    -------
        (M - Mhat)' * W * (M - Mhat): weighted distance between theoretical and empirical moments
    """

    # compute theoretical moments, then compute distance from empirical moments
    M = theoretical_moments(params, data)
    return (M - Mhat)' * W * (M - Mhat)
end

# First stage: minimize objective using equal weights (i.e. W = Identity matrix)
Mhat, T = empirical_moments(reslogwage)
W = I(length(Mhat))
r = optimize(params->GMM_objective(params, Mhat, reslogwage, W), [0.8, 0.1, 0.1, 0.1], LBFGS())

# Second stage: use first-stage GMM estimates to compute W according to optimal weighting matrix formula.
# To avoid small sample bias, extract diagonal to minimize Diagonally Weighted Minimum Distance in the next
# iteration of GMM (following Blundell, Pistaferri, and Preston (2008)).
M = theoretical_moments(r.minimizer, reslogwage)
W = I(length(Mhat)) .* diag(pinv((1 / T) * M * M'))
r = optimize(params->GMM_objective(params, Mhat, reslogwage, W), [0.8, 0.1, 0.1, 0.1], LBFGS())
rho, var_alpha, var_nu, var_eta = r.minimizer

# Save estimates to access later
const rho_GMM, var_alpha_GMM, var_nu_GMM, var_eta_GMM = rho, var_alpha, var_nu, var_eta

end