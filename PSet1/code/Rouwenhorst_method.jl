# Approximate income with a Finite-State Markov Process, discretizing using the Rouwenhorst method.

module Rouwenhorst_method
export Rouwenhorst

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

end