"""
primal and dual sub-problems
"""

function post_dc_primal(data::Dict{String,Any}, scenarios, model=Model())   
    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]
    
    gen_keys = collect(keys(ref[:gen]))
    br_keys = collect(keys(ref[:branch]))
    
    for i in 1:length(gen_keys)
        gen_id = gen_keys[i]
        if scenarios[i] == 1
            # delete!(data["gen"], string(gen_id))
            data["gen"][string(gen_id)]["gen_status"] = 0
        end
    end

    for i in 1:length(br_keys)
        branch_id = br_keys[i]
        if scenarios[i+length(gen_keys)] == 1
            # delete!(data["branch"], string(branch_id))
            data["branch"][string(branch_id)]["br_status"] = 0
        end
    end

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
    tmin_cons = JuMP.ConstraintRef[]
    tmax_cons = JuMP.ConstraintRef[]
    dclb_cons = JuMP.ConstraintRef[]
    dcub_cons = JuMP.ConstraintRef[]
    vamin_cons = JuMP.ConstraintRef[]
    vamax_cons = JuMP.ConstraintRef[]
    loadshed_cons = JuMP.ConstraintRef[]

    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        bus_gens = ref[:bus_gens][i]

        push!(kcl_cons, @constraint(model, sum(p_expr[a] for a in bus_arcs) - sum(pg[g] for g in bus_gens)  - ld[i]*bus["pd"] == - bus["pd"]))
    end

    for (i, gen) in ref[:gen]
        push!(pgmin_cons, @constraint(model, pg[i] >= 0))
        push!(pgmax_cons, @constraint(model, pg[i] <= gen["pmax"]))
    end

    for (l,i,j) in ref[:arcs_from]
        push!(tmin_cons, @constraint(model, p[(l,i,j)] >= -ref[:branch][l]["rate_a"]))
        push!(tmax_cons, @constraint(model, p[(l,i,j)] <= ref[:branch][l]["rate_a"]))
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PMs.calc_branch_y(branch)

        push!(dclb_cons, @constraint(model, p_fr + b*(va_fr - va_to) >= 0))
        push!(dcub_cons, @constraint(model, p_fr + b*(va_fr - va_to) <= 0))
    end

    for (i, branch) in ref[:branch]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        push!(vamin_cons, @constraint(model, va_fr - va_to >= branch["angmin"]))
        push!(vamax_cons, @constraint(model, va_fr - va_to <= branch["angmax"]))
    end

    for (i, bus) in ref[:bus]
        push!(loadshed_cons, @constraint(model, ld[i] <= 1))
    end


    @objective(model, Min, sum(ref[:bus][i]["pd"] * ld[i] for i in keys(ref[:bus])))
    
    return model
    
end

function post_dc_dual(data::Dict{String,Any}, scenarios, model=Model())
    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]
    
    gen_keys = collect(keys(ref[:gen]))
    br_keys = collect(keys(ref[:branch]))

    for i in 1:length(gen_keys)
        gen_id = gen_keys[i]
        if scenarios[i] == 1
            # delete!(data["gen"], string(gen_id))
            data["gen"][string(gen_id)]["gen_status"] = 0
        end
    end

    for i in 1:length(br_keys)
        branch_id = br_keys[i]
        if scenarios[i+length(gen_keys)] == 1
            # delete!(data["branch"], string(branch_id))
            data["branch"][string(branch_id)]["br_status"] = 0
        end
    end

    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]

    @variable(model, kcl[i in keys(ref[:bus])])
    @variable(model, pgmin[i in keys(ref[:gen])] >= 0)
    @variable(model, pgmax[i in keys(ref[:gen])] >= 0)
    @variable(model, tmin[i in keys(ref[:branch])] >= 0)
    @variable(model, tmax[i in keys(ref[:branch])] >= 0)
    @variable(model, dclb[i in keys(ref[:branch])] >= 0)
    @variable(model, dcub[i in keys(ref[:branch])] >= 0)
    @variable(model, vamin[i in keys(ref[:branch])] >= 0)
    @variable(model, vamax[i in keys(ref[:branch])] >= 0)
    @variable(model, loadshed[i in keys(ref[:bus])] >= 0)

    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))
    
    va_cons = JuMP.ConstraintRef[]
    pg_cons = JuMP.ConstraintRef[]
    p_cons = JuMP.ConstraintRef[]
    ld_cons = JuMP.ConstraintRef[]


    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i] 
        push!(va_cons, 
              @constraint(model, 
                          sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dclb[l] - dcub[l]) + p_mag[(l,f,t)] * (vamin[l] - vamax[l]) for (l,f,t) in bus_arcs) == 0))
    end

    for (i, gen) in ref[:gen]
        push!(pg_cons, @constraint(model, -kcl[gen["gen_bus"]] + pgmin[i] - pgmax[i] == 0))
    end

    for (l,i,j) in ref[:arcs_from]
        push!(p_cons, @constraint(model, kcl[i] - kcl[j] + tmin[l] - tmax[l] + dclb[l] - dcub[l] == 0))
    end

    for (i, bus) in ref[:bus]
        push!(ld_cons, @constraint(model, -bus["pd"]*kcl[i] - loadshed[i] <= bus["pd"]))
    end
    

    @objective(model, Max, 
               sum( -ref[:bus][i]["pd"] * kcl[i] for i in keys(ref[:bus]) ) +
               sum( -ref[:gen][i]["pmax"] * pgmax[i] for i in keys(ref[:gen]) ) +
               sum( -ref[:branch][l]["rate_a"] * tmin[l] - ref[:branch][l]["rate_a"] * tmax[l] for l in keys(ref[:branch]) ) +
               sum( 0 + 0 for i in keys(ref[:branch]) ) + 
               sum( ref[:branch][l]["angmin"] * vamin[l] - ref[:branch][l]["angmax"] * vamax[l] for l in keys(ref[:branch]) ) +
               sum( -loadshed[i] for i in keys(ref[:bus]) )
              )

    return model

end

function post_dc_kkt(data::Dict{String,Any}, scenarios, model=Model())
    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]
    
    gen_keys = collect(keys(ref[:gen]))
    br_keys = collect(keys(ref[:branch]))
    
    for i in 1:length(gen_keys)
        gen_id = gen_keys[i]
        if scenarios[i] == 1
            # delete!(data["gen"], string(gen_id))
            data["gen"][string(gen_id)]["gen_status"] = 0
        end
    end

    for i in 1:length(br_keys)
        branch_id = br_keys[i]
        if scenarios[i+length(gen_keys)] == 1
            # delete!(data["branch"], string(branch_id))
            data["branch"][string(branch_id)]["br_status"] = 0
        end
    end

    ref = PMs.build_ref(data)
    
    ref = ref[:nw][0]

    # primal variables
    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, pg[i in keys(ref[:gen])])
    @variable(model, p[(l,i,j) in ref[:arcs_from]])
    @variable(model, ld[i in keys(ref[:bus])] >= 0)

    # dual variables
    @variable(model, kcl[i in keys(ref[:bus])])
    @variable(model, pgmin[i in keys(ref[:gen])] >= 0)
    @variable(model, pgmax[i in keys(ref[:gen])] >= 0)
    @variable(model, tmin[i in keys(ref[:branch])] >= 0)
    @variable(model, tmax[i in keys(ref[:branch])] >= 0)
    @variable(model, dclb[i in keys(ref[:branch])] >= 0)
    @variable(model, dcub[i in keys(ref[:branch])] >= 0)
    @variable(model, vamin[i in keys(ref[:branch])] >= 0)
    @variable(model, vamax[i in keys(ref[:branch])] >= 0)
    @variable(model, loadshed[i in keys(ref[:bus])] >= 0)

    # auxiliary definitions for formulating the constraints
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))

    # primal constraint references
    kcl_cons = JuMP.ConstraintRef[]
    pgmin_cons = JuMP.ConstraintRef[]
    pgmax_cons = JuMP.ConstraintRef[]
    tmin_cons = JuMP.ConstraintRef[]
    tmax_cons = JuMP.ConstraintRef[]
    dclb_cons = JuMP.ConstraintRef[]
    dcub_cons = JuMP.ConstraintRef[]
    vamin_cons = JuMP.ConstraintRef[]
    vamax_cons = JuMP.ConstraintRef[]
    loadshed_cons = JuMP.ConstraintRef[]

    # dual constraint references
    va_cons = JuMP.ConstraintRef[]
    pg_cons = JuMP.ConstraintRef[]
    p_cons = JuMP.ConstraintRef[]
    ld_cons = JuMP.ConstraintRef[]

    # other constraint references
    pdeq_cons = JuMP.ConstraintRef[]
    
    # useful expressions
    primalobj_expr = Any
    dualobj_expr = Any

    # primal constraints
    # (a) kcl
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        bus_gens = ref[:bus_gens][i]

        push!(kcl_cons, @constraint(model, sum(p_expr[a] for a in bus_arcs) - sum(pg[g] for g in bus_gens)  - ld[i]*bus["pd"] == - bus["pd"]))
    end

    # (b) generation limits
    for (i, gen) in ref[:gen]
        push!(pgmin_cons, @constraint(model, pg[i] >= 0))
        push!(pgmax_cons, @constraint(model, pg[i] <= gen["pmax"]))
    end

    # (c) thermal limits
    for (l,i,j) in ref[:arcs_from]
        push!(tmin_cons, @constraint(model, p[(l,i,j)] >= -ref[:branch][l]["rate_a"]))
        push!(tmax_cons, @constraint(model, p[(l,i,j)] <= ref[:branch][l]["rate_a"]))
    end

    # (d) dc power flow
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])

        p_fr = p[f_idx]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        g, b = PMs.calc_branch_y(branch)

        push!(dclb_cons, @constraint(model, p_fr + b*(va_fr - va_to) >= 0))
        push!(dcub_cons, @constraint(model, p_fr + b*(va_fr - va_to) <= 0))
    end

    # (e) va bounds
    for (i, branch) in ref[:branch]
        va_fr = va[branch["f_bus"]]
        va_to = va[branch["t_bus"]]

        push!(vamin_cons, @constraint(model, va_fr - va_to >= branch["angmin"]))
        push!(vamax_cons, @constraint(model, va_fr - va_to <= branch["angmax"]))
    end

    # (f) load shedding limit
    for (i, bus) in ref[:bus]
        push!(loadshed_cons, @constraint(model, ld[i] <= 1))
    end
    
    # dual constraints
    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i] 
        push!(va_cons, 
              @constraint(model, 
                          sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dclb[l] - dcub[l]) + p_mag[(l,f,t)] * (vamin[l] - vamax[l]) for (l,f,t) in bus_arcs) == 0))
    end
    
    # (b) constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
        push!(pg_cons, @constraint(model, -kcl[gen["gen_bus"]] + pgmin[i] - pgmax[i] == 0))
    end

    # (c) constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        push!(p_cons, @constraint(model, kcl[i] - kcl[j] + tmin[l] - tmax[l] + dclb[l] - dcub[l] == 0))
    end

    # (d) constraints corresponding to primal variable ld
    for (i, bus) in ref[:bus]
        push!(ld_cons, @constraint(model, -bus["pd"]*kcl[i] - loadshed[i] <= bus["pd"]))
    end


    # primal objective == dual 
    @expression(model, primalobj_expr, sum( ref[:bus][i]["pd"] * ld[i] for i in keys(ref[:bus]) ) )
    @expression(model, dualobj_expr, 
                                sum( -ref[:bus][i]["pd"] * kcl[i] for i in keys(ref[:bus]) ) +
                                sum( -ref[:gen][i]["pmax"] * pgmax[i] for i in keys(ref[:gen]) ) +
                                sum( -ref[:branch][l]["rate_a"] * tmin[l] - ref[:branch][l]["rate_a"] * tmax[l] for l in keys(ref[:branch]) ) +
                                sum( 0 + 0 for i in keys(ref[:branch]) ) + 
                                sum( ref[:branch][l]["angmin"] * vamin[l] - ref[:branch][l]["angmax"] * vamax[l] for l in keys(ref[:branch]) ) +
                                sum( -loadshed[i] for i in keys(ref[:bus]) )
                               )
    
    model.ext[:primalobj_expr] = primalobj_expr

    push!(pdeq_cons, @constraint(model, primalobj_expr - dualobj_expr == 0))

    return model

end

