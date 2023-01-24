using .GMM, .Rouwenhorst_method, Interpolations, Plots, LaTeXStrings

default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

N = 5
y_grid, Pi = Rouwenhorst(rho, sqrt(var_eta), N)
# y_grid, Theta = Rouwenhorst(0.985, sqrt(0.0346), N) # params from top 1% paper

# prob_mat = [0.969909 0.029317 0.000332 0.000002 0.000000
#             0.007329 0.970075 0.021989 0.000166 0.000000
#             0.000055 0.014659 0.970130 0.014659 0.000055
#             0.000000 0.000166 0.021989 0.970075 0.007329
#             0.000000 0.000002 0.000332 0.029317 0.969909]

# row_sum = sum(prob_mat, dims=2)
# prob_mat = prob_mat ./ row_sum

beta = 0.96
r = 0.02
gamma = 2
a_min = 0

u(c) = c^(1 - gamma) / (1 - gamma)
u_prime(c) = c^(-gamma)
c_hat(B) = B.^(-1 / gamma)

G_a  = [2500 * ((1 + 0.08)^(i - 1) - 1) / ((1 + 0.08)^99 - 1) for i in 1:100] # [2500 * ((1 + 0.08)^(i - 1) - 1) / ((1 + 0.08)^99 - 1) for i in 1:100]
G_y = exp.(y_grid)
c = r * G_a .+ G_y'
tol = 1e-9
dist = 1000

while dist > tol
    B = zeros(length(G_a), length(G_y))
    for (i, ap) in enumerate(G_a)
        for (j, y) in enumerate(G_y)
            B[i, j] = beta * (1 + r) * sum(Pi[j, :] .* u_prime.(c[i, :]))
        end
    end

    a_star = (c_hat(B) .+ G_a .- G_y') / (1 + r)
    a_star_min = a_star[1, :]

    c_old = copy(c)
    for (j, y) in enumerate(G_y)
        interp = LinearInterpolation(a_star[:, j], c_hat(B)[:, j], extrapolation_bc=Line())
        c[:, j] = interp.(G_a)

        c[G_a .< a_star_min[j], j] = (1 + r) .* G_a[G_a .< a_star_min[j]] .+ y .- a_min
    end
    dist = maximum(abs.(c .- c_old))
end

plot(G_a, c[:, 3])
plot!(G_a, c[:, 1])
plot!(G_a, c[:, 5])



plot(G_a[G_a .< 1], c[G_a .< 1, 3])
plot!(G_a[G_a .< 1], c[G_a .< 1, 1])
plot!(G_a[G_a .< 1], c[G_a .< 1, 5])

# plot(x, f.(x), 
#     label=L"x^2",
#     xlabel=L"x",
#     ylabel=L"x^2",
#     ylims = (0,25))
# savefig("intro/x-squared-jl.pdf")
