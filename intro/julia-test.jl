using Plots
using LaTeXStrings
using Optim

default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

f(x) = x^2;
x = LinRange(-5, 5, 100);

plot(x, f.(x), 
    label=L"x^2",
    xlabel=L"x",
    ylabel=L"x^2",
    ylims = (0,25))
savefig("intro/x-squared-jl.pdf")

res = optimize(f, -5.0, 5.0)
println("Minimizer: ",  round(res.minimizer, digits=3, base=10));
println("Minimum: ",  round(res.minimum, digits=3, base=10));
