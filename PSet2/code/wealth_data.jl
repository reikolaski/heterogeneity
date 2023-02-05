using DataFrames, CSV, Plots, LaTeXStrings, StatsBase
default(fontfamily="Computer Modern", framestyle=:box, label=nothing, grid=false, legendfontsize=10)
scalefontsizes(1)

data = hcat(DataFrame(CSV.File("data/SCFP2019.csv")), DataFrame(CSV.File("data/people.csv")))
# WGT, NETWORTH, KIDS

# test = data[:, [:WGT, :KIDS, :FAMSTRUCT, :NETWORTH, :PEOPLE, :PEOPLE_HHL]]
# test.ADULTS1 .= ifelse.(test[:, :FAMSTRUCT] .<= 3, 1, 2)
# test.ADULTS2 .= test[:, :PEOPLE] .- test[:, :KIDS]
# test.ADULTS3 .= test[:, :PEOPLE_HHL] .- test[:, :KIDS]
# describe(test)

scf = data[:, [:WGT, :AGE, :AGECL, :FAMSTRUCT, :KIDS, :NETWORTH, :PEOPLE]]
scf.ADULTS .= scf[:, :PEOPLE] .- scf[:, :KIDS]
scf.ADULTS .= ifelse.(scf.ADULTS .< 1, 1, scf.ADULTS)
# scf.ADULTS .= ifelse.(scf[:, :FAMSTRUCT] .<= 3, 1, 2)
scf.EFFECTIVE_PEOPLE .= 1 .+ 0.7(scf.ADULTS .- 1) .+ 0.5(scf.KIDS)
scf.EFFECTIVE_NETWORTH .= scf.NETWORTH ./ scf.EFFECTIVE_PEOPLE

scf.AGE_BIN = scf.AGECL .+ 1
scf.AGE_BIN .= ifelse.(scf.AGE .< 25, 1, scf.AGE_BIN)

avg_lc_networth = combine(groupby(scf, :AGE_BIN), [:EFFECTIVE_NETWORTH, :WGT] => ((x, y) -> mean(x, fweights(y))) => :weighted_mean)
avg_lc_networth.BIN_DEF = ["<25", "25-34", "35-44", "45-54", "55-64", "65-74", "75+"]
plot(avg_lc_networth.AGE_BIN, avg_lc_networth.weighted_mean / 1000,
     xticks = (1:7, avg_lc_networth.BIN_DEF),
     xlabel="Age Cohort",
     ylabel="Average Net Worth (Thousands)")
savefig("plots/average_cohort_net_worth.pdf")

# sum(scf.EFFECTIVE_NETWORTH[scf.AGE .== 60] .* scf.WGT[scf.AGE .== 60]) / sum(scf.WGT[scf.AGE .== 60])
# sum(scf.EFFECTIVE_NETWORTH[scf.AGE_BIN .== 5] .* scf.WGT[scf.AGE_BIN .== 5]) / sum(scf.WGT[scf.AGE_BIN .== 5])
scf60 = scf[scf.AGE .== 60, :]
sort!(scf60, [:EFFECTIVE_NETWORTH])
scf60.CUMWGT = cumsum(scf60.WGT)
scf60.EFFECTIVE_NETWORTH[findfirst(x -> x > (maximum(scf60.CUMWGT) / 2), scf60.CUMWGT)]


# sort!(scf, [:EFFECTIVE_NETWORTH])
# scf.CUMWGT = cumsum(scf.WGT)
# scf.EFFECTIVE_NETWORTH[findfirst(x -> x > (maximum(scf.CUMWGT) / 2), scf.CUMWGT)]
