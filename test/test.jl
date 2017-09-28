using PowerModels
using JuMP
using CPLEX
include("../src/subproblem.jl")

PMs = PowerModels

function gen_model(;formulation="primal", case=3)
    (case == 3) ? (data = PMs.parse_file("../data/pglib_opf_case3_lmbd_nk.m")) : (data = PMs.parse_file("../data/pglib_opf_case24_ieee_rts.m"))

    if formulation == "primal"
        m = post_dc_primal(data, Model(solver=CplexSolver()))
        return m
    end

    if formulation == "dual"
        m = post_dc_dual(data, Model(solver=CplexSolver()))
        return m
    end

    if formulation == "kkt"
        m = post_dc_kkt(data, Model(solver=CplexSolver()))
        return m
    end
end
