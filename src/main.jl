using JuMP
using PowerModels
using CPLEX

PMs = PowerModels

include("subproblem.jl")

data = PMs.parse_file("../data/pglib_opf_case24_ieee_rts.m")
model = post_dc_primal(data, Model(solver=CplexSolver()))

status = solve(model)

model = post_dc_dual(data, Model(solver=CplexSolver()))

status = solve(model)
