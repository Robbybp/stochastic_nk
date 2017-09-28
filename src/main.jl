using JuMP
using PowerModels
using CPLEX

PMs = PowerModels

include("subproblem.jl")

data = PMs.parse_file("../data/pglib_opf_case24_ieee_rts.m")
m_primal = post_dc_primal(data, Model(solver=CplexSolver()))

status = solve(m_primal)

m_dual = post_dc_dual(data, Model(solver=CplexSolver()))

status = solve(m_dual)

@assert getobjectivevalue(m_primal) == getobjectivevalue(m_dual)
