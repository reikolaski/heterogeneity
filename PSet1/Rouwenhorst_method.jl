module Rouwenhorst_method
export Rouwenhorst

function Rouwenhorst(rho, sigma, N)
    p = (1 + rho) / 2
    q = p
    Psi = sqrt(N - 1) * sigma

    Theta_m1 = [[p 1-p];
                [1-q q]]
    for n in 3:N
        global Theta = p .* [[Theta_m1 zeros(n-1, 1)]; [zeros(1, n-1) 0]] +
                        (1 - p) .* [[zeros(n-1, 1) Theta_m1]; [0 zeros(1, n-1)]] +
                        (1 - q) .* [[zeros(1, n-1) 0]; [Theta_m1 zeros(n-1, 1)]] +
                        q .* [[0 zeros(1, n-1)]; [zeros(n-1, 1) Theta_m1]]
        Theta[2:end-1, :] /= 2
        Theta_m1 = Theta
    end
    
    y_grid = LinRange(-Psi, Psi, N)
    return y_grid, Theta
end

end