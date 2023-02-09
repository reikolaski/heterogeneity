module SCF_wealth_moments
export SCF_mean_wealth_age_60, SCF_median_wealth_age_60

using DataFrames, CSV, Plots, LaTeXStrings, StatsBase
default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

# Load in SCF data
data = hcat(DataFrame(CSV.File("data/SCFP2019.csv")), DataFrame(CSV.File("data/people.csv")))

# Compute effective net worth of household head using Oxford equivalence scale
scf = data[:, [:WGT, :AGE, :AGECL, :FAMSTRUCT, :KIDS, :NETWORTH, :PEOPLE]]
scf.ADULTS .= scf[:, :PEOPLE] .- scf[:, :KIDS]
scf.ADULTS .= ifelse.(scf.ADULTS .< 1, 1, scf.ADULTS)
scf.EFFECTIVE_PEOPLE .= 1 .+ 0.7(scf.ADULTS .- 1) .+ 0.5(scf.KIDS)
scf.EFFECTIVE_NETWORTH .= scf.NETWORTH ./ scf.EFFECTIVE_PEOPLE

# Create age bins
scf.AGE_BIN = scf.AGECL .+ 1
scf.AGE_BIN .= ifelse.(scf.AGE .< 25, 1, scf.AGE_BIN)

# Compute weighted average net worth within each age bin
avg_lc_networth = combine(groupby(scf, :AGE_BIN), [:EFFECTIVE_NETWORTH, :WGT] => ((x, y) -> mean(x, fweights(y))) => :weighted_mean)
avg_lc_networth.BIN_DEF = ["<25", "25-34", "35-44", "45-54", "55-64", "65-74", "75+"]
plot(avg_lc_networth.AGE_BIN, avg_lc_networth.weighted_mean / 1000,
     xticks = (1:7, avg_lc_networth.BIN_DEF),
     xlabel="Age Cohort",
     ylabel="Average Net Worth (Thousands)")
savefig("plots/average_cohort_net_worth.pdf")

# Compute average net worth at age 60
const SCF_mean_wealth_age_60 = sum(scf.EFFECTIVE_NETWORTH[scf.AGE .== 60] .* scf.WGT[scf.AGE .== 60]) / sum(scf.WGT[scf.AGE .== 60])
println("\nMean net worth at age 60: $SCF_mean_wealth_age_60")

# Compute average net worth at age 60
scf60 = scf[scf.AGE .== 60, :]
sort!(scf60, [:EFFECTIVE_NETWORTH])
scf60.CUMWGT = cumsum(scf60.WGT)
const SCF_median_wealth_age_60 = scf60.EFFECTIVE_NETWORTH[findfirst(x -> x > (maximum(scf60.CUMWGT) / 2), scf60.CUMWGT)]
println("\nMedian net worth at age 60: $SCF_median_wealth_age_60")

end