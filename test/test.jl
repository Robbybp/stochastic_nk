using PowerModels
using JuMP
using CPLEX
include("../src/subproblem.jl")

PMs = PowerModels

function gen_model(;primal=true)
    data = PMs.parse_file("../data/pglib_opf_case24_ieee_rts.m")
    (primal) ? (m = post_dc_primal(data, Model(solver=CplexSolver()))) : (m = post_dc_dual(data, Model(solver=CplexSolver())))
    return m
end
