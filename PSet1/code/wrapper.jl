# Runs all Julia code for Econ 237 Problem Set 1 (after cleaning data in Stata)
# Working directory should be set as "/PSet1/"

# Load all packages
using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("StatFiles")
Pkg.add("Optim")
Pkg.add("Interpolations")
Pkg.add("Plots")
Pkg.add("LaTeXStrings")

# Run programs
include("GMM.jl")
include("Rouwenhorst_method.jl")
include("IFP.jl")