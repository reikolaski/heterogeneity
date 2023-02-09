# Runs all Julia code for Econ 237 Problem Set 2

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
Pkg.add("Distributions")
Pkg.add("DataStructures")

cd("/Users/reiko/Documents/heterogeneity/PSet2/")

# Run programs
include("wealth_data.jl")
include("EGM_functions.jl")
include("EGM.jl")