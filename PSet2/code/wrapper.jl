# Runs all Julia code for Econ 237 Problem Set 2
# Working directory should be set as "/PSet2/"

# Load all packages
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("StatFiles")
Pkg.add("Optim")
Pkg.add("Interpolations")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")
Pkg.add("StatsBase")

cd("/Users/reiko/Documents/heterogeneity/PSet2/")


# Run programs
# include("GMM.jl")
include("Rouwenhorst_method.jl")
# include("IFP.jl")