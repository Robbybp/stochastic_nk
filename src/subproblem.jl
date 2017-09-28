"""
primal and dual sub-problems
"""

function post_dc_primal(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]

    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, p[(l,i,j) in ref[:arcs_from]])
    @variable(model, ld[i in keys(ref[:bus])] >= 0)

    
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
    
    kcl_cons = JuMP.ConstraintRef[]
    pgmin_cons = JuMP.ConstraintRef[]
    pgmax_cons = JuMP.ConstraintRef[]
    thermal_min_cons = JuMP.ConstraintRef[]
    thermal_max_cons = JuMP.ConstraintRef[]
    dc_lb_cons = JuMP.ConstraintRef[]
    dc_ub_cons = JuMP.ConstraintRef[]
    load_shed_cons = JuMP.ConstraintRef[]

    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        bus_gens = ref[:bus_gens][i]

        push!(kcl_cons, @constraint(model, sum(p_expr[a] for a in bus_arcs) - sum(pg[g] for g in bus_gens)  - ld[i]*bus["pd"] == - bus["pd"]))
    end

    for (i, gen) in ref[:gen]
        push!(pgmin_cons, @constraint(model, pg[i] >= gen["pmin"]))
        push!(pgmax_cons, @constraint(model, pg[i] <= gen["pmax"]))
    end

    for (l,i,j) in ref[:arcs_from]
        push!(thermal_min_cons, @constraint(model, p[(l,i,j)] >= -ref[:branch][l]["rate_a"]))
        push!(thermal_max_cons, @constraint(model, p[(l,i,j)] <= ref[:branch][l]["rate_a"]))
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PMs.calc_branch_y(branch)

        push!(dc_lb_cons, @constraint(model, p_fr + b*(va_fr - va_to) >= 0))
        push!(dc_ub_cons, @constraint(model, p_fr + b*(va_fr - va_to) <= 0))
    end

    for (i, bus) in ref[:bus]
        push!(load_shed_cons, @constraint(model, ld[i] <= 1))
    end


    @objective(model, Min, sum(ref[:bus][i]["pd"] * ld[i] for i in keys(ref[:bus])))
    
    return model
    
end

function post_dc_dual(data::Dict{String,Any}, model=Model())
    ref = PMs.build_ref(data)

    ref = ref[:nw][0]

    @variable(model, kcl_dual[i in keys(ref[:bus])])
    @variable(model, pgmin_dual[i in keys(ref[:gen])] >= 0)
    @variable(model, pgmax_dual[i in keys(ref[:gen])] >= 0)
    @variable(model, thermal_min_dual[i in keys(ref[:branch])] >= 0)
    @variable(model, thermal_max_dual[i in keys(ref[:branch])] >= 0)
    @variable(model, dc_lb_dual[i in keys(ref[:branch])] >= 0)
    @variable(model, dc_ub_dual[i in keys(ref[:branch])] >= 0)
    @variable(model, load_shed_dual[i in keys(ref[:bus])] >= 0)

    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))
    
    va_dual_cons = JuMP.ConstraintRef[]
    pg_dual_cons = JuMP.ConstraintRef[]
    p_dual_cons = JuMP.ConstraintRef[]
    ld_dual_cons = JuMP.ConstraintRef[]


    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i] 
        push!(va_dual_cons, 
              @constraint(model, 
              sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dc_lb_dual[l] + dc_ub_dual[l]) for (l,f,t) in bus_arcs) == 0))
    end

    for (i, gen) in ref[:gen]
        push!(pg_dual_cons, @constraint(model, -kcl_dual[gen["gen_bus"]] + pgmin_dual[i] - pgmax_dual[i] == 0))
    end

    for (l,i,j) in ref[:arcs_from]
        push!(p_dual_cons, @constraint(model, kcl_dual[i] - kcl_dual[j] + thermal_min_dual[l] - thermal_max_dual[l] + dc_lb_dual[l] - dc_ub_dual[l] == 0))
    end

    for (i, bus) in ref[:bus]
        push!(ld_dual_cons, @constraint(model, -bus["pd"]*kcl_dual[i] - load_shed_dual[i] <= bus["pd"]))
    end
    

    @objective(model, Max, 
               sum( -ref[:bus][i]["pd"] * kcl_dual[i] for i in keys(ref[:bus]) ) +
               sum( ref[:gen][i]["pmin"] * pgmin_dual[i] - ref[:gen][i]["pmax"] * pgmax_dual[i] for i in keys(ref[:gen]) ) +
               sum( -ref[:branch][l]["rate_a"] * thermal_min_dual[l] - ref[:branch][l]["rate_a"] * thermal_max_dual[l] for l in keys(ref[:branch]) ) +
               sum( 0 + 0 for i in keys(ref[:branch]) ) + 
               sum( -load_shed_dual[i] for i in keys(ref[:bus]) )
              )

    return model

end

