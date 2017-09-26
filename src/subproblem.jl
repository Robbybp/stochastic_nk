"""
primal and dual sub-problems
"""

function post_dc_primal(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, ref[:gen][i]["pmin"] <= pg[i in keys(ref[:gen])] <= ref[:gen][i]["pmax"])
    @variable(model, -ref[:branch][l]["rate_a"] <= p[(l,i,j) in ref[:arcs_from]] <= ref[:branch][l]["rate_a"])

    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))

    ld_max = Dict(i => 0.0 for i in keys(ref[:bus]))
    for i in keys(ref[:bus])
    ld_max[i] = (ref[:bus][i]["pd"] > 0) ? 1.0 : 0.0
    end

    @variable(model, 0 <= ld[i in keys(ref[:bus])] <= ld_max[i])

    @objective(model, Min, sum(ref[:bus][i]["pd"] * ld[i] for i in keys(ref[:bus])))
    
    for (i, bus) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        bus_gens = ref[:bus_gens][i]

        @constraint(model, sum(p_expr[a] for a in bus_arcs) == sum(pg[g] for g in bus_gens) - bus["pd"] + ld[i]*bus["pd"])
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PMs.calc_branch_y(branch)

        @constraint(model, p_fr == -b*(va_fr - va_to))
    end
    return model
    
end

function post_dc_dual(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)

    ref = ref[:nw][0]

    # dual variable corresponding to kcl
    @variable(model, α[i in keys(ref[:branch])])

    # dual variable corresponding to pgmax
    #@variable(model, β⁺[i in 

    return model

end

