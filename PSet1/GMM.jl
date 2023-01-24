module GMM
export rho, var_alpha, var_nu, var_eta

using CSV, DataFrames, StatFiles, Statistics, LinearAlgebra, Optim

reslogwage = DataFrame(CSV.File("residual_logwage.csv"))

function empirical_moments(data)
    T = ncol(data)
    cov_mat = zeros(T, T)
    for (i, c1) in enumerate(eachcol(data))
        for (j, c2) in enumerate(eachcol(data))
            sx, sy = skipmissings(c1, c2)
            cov_mat[i, j] = cov(collect(sx), collect(sy))
        end
    end
    return cov_mat[triu!(trues(size(cov_mat)))], T
end

function theoretical_moments(params, data)
    rho, var_alpha, var_nu, var_eta = params
    T = ncol(data)
    cov_mat = zeros(T, T)
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
    return cov_mat[triu!(trues(size(cov_mat)))]
end

function GMM_objective(params, Mhat, reslogwage, W)
    M = theoretical_moments(params, reslogwage)
    return (M - Mhat)' * W * (M - Mhat)
end

Mhat, T = empirical_moments(reslogwage)
W = I(length(Mhat))
r = optimize(params->GMM_objective(params, Mhat, reslogwage, W), [0.8, 0.1, 0.1, 0.1], LBFGS())

M = theoretical_moments(r.minimizer, reslogwage)
W = I(length(Mhat)) .* diag(pinv((1 / T) * M * M'))
r = optimize(params->GMM_objective(params, Mhat, reslogwage, W), [0.8, 0.1, 0.1, 0.1], LBFGS())
const rho, var_alpha, var_nu, var_eta = r.minimizer

end