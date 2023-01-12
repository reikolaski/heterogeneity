# using Plots, Random

f(x) = x + 1

function plot_results()
    x = 1:5
    y = f.(x)
    plot(x, y)
    print(y)
end

plot_results()

using LinearAlgebra, Statistics, Plots, LaTeXStrings
n = 1000000
eps = zeros(n)
for i in eachindex(eps)
    eps[i] = randn()
end

f(x) = x^2
generatedata(n) = f.(randn(n))
data = generatedata(5)