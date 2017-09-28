
function create_full_model(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}; model=Model())
    
    ref = ref[:nw][0]
    numscenarios = config["numscenarios"]

    # interdiction variables
    @variable(model, x[i in keys(ref[:branch])], Bin)
    @variable(model, y[i in keys(ref[:gen])], Bin)

    # primal variables
    @variable(model, va[i in keys(ref[:bus]), s in 1:numscenarios])
    @variable(model, pg[i in keys(ref[:gen]), s in 1:numscenarios])
    @variable(model, p[(l,i,j) in ref[:arcs_from], s in 1:numscenarios])
    @variable(model, ld[i in keys(ref[:bus]), s in 1:numscenarios] >= 0)

    # dual variables
    @variable(model, kcl[i in keys(ref[:bus]), s in 1:numscenarios])
    @variable(model, pgmin[i in keys(ref[:gen]), s in 1:numscenarios] >= 0)
    @variable(model, pgmax[i in keys(ref[:gen]), s in 1:numscenarios] >= 0)
    @variable(model, tmin[i in keys(ref[:branch]), s in 1:numscenarios] >= 0)
    @variable(model, tmax[i in keys(ref[:branch]), s in 1:numscenarios] >= 0)
    @variable(model, dclb[i in keys(ref[:branch]), s in 1:numscenarios] >= 0)
    @variable(model, dcub[i in keys(ref[:branch]), s in 1:numscenarios] >= 0)
    @variable(model, loadshed[i in keys(ref[:bus]), s in 1:numscenarios] >= 0)

    # auxiliary definitions for formulating the constraints
    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))

    # useful expressions
    primalobj_expr = Any[]
    dualobj_expr = Any[]

    # master constraints
    @constraint(model, sum(x) + sum(y) == config["budget"])
    (config["interdict"] == "l") && @constraint(model, sum(y) == 0)

    # primal constraints
    # (a) kcl
    for (b, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][b]
        bus_gens = ref[:bus_gens][b]
        
        for s in 1:numscenarios
            p_expr = Dict([((l,i,j), 1.0*p[(l,i,j),s]) for (l,i,j) in ref[:arcs_from]])
            p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j),s]) for (l,i,j) in ref[:arcs_from]]))
            @constraint(model, sum(p_expr[a] for a in bus_arcs) - sum(pg[g,s] for g in bus_gens)  - ld[b,s]*bus["pd"] == - bus["pd"])
        end
        
    end

    # (b) generation limits
    for i in 1:length(keys(ref[:gen]))
        gen_id = collect(keys(ref[:gen]))[i]
        gen = ref[:gen][gen_id]
        @constraint(model, [s=1:numscenarios], pg[gen_id,s] >= (1-y[gen_id]*scenarios[s,i]) * gen["pmin"])
        @constraint(model, [s=1:numscenarios], pg[gen_id,s] <= (1-y[gen_id]*scenarios[s,i]) * gen["pmax"])
    end
    
    # (c) thermal limits
    for k in 1:length(ref[:arcs_from])
        l, i, j = ref[:arcs_from][k]
        col_index = length(keys(ref[:gen])) + k
        @constraint(model, [s=1:numscenarios], p[(l,i,j),s] >= -ref[:branch][l]["rate_a"] * (1-x[l]*scenarios[s,col_index]))
        @constraint(model, [s=1:numscenarios], p[(l,i,j),s] <= ref[:branch][l]["rate_a"] * (1-x[l]*scenarios[s,col_index]))
    end
    
    # dc power flow
    for k in 1:length(keys(ref[:branch]))
        i = collect(keys(ref[:branch]))[k]
        branch = ref[:branch][i]
        col_index = length(keys(ref[:gen])) + k

        f_idx = (i, branch["f_bus"], branch["t_bus"])

        g, b = PMs.calc_branch_y(branch)
    
        for s in 1:numscenarios
            p_fr = p[f_idx,s]
            va_fr = va[branch["f_bus"],s]
            va_to = va[branch["t_bus"],s]

            @constraint(model, p_fr + b*(va_fr - va_to) >= ref[:off_angmin]*x[i]*scenarios[s,col_index])
            @constraint(model, p_fr + b*(va_fr - va_to) <= ref[:off_angmax]*x[i]*scenarios[s,col_index])
        end
    end

    # load shedding limit
    for (i, bus) in ref[:bus]
        @constraint(model, [s=1:numscenarios], ld[i,s] <= 1)
    end

    println(model)
    quit()
    
    # dual constraints
    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i] 
        push!(va_cons, 
              @constraint(model, 
              sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (tmin[l] + tmax[l]) for (l,f,t) in bus_arcs) == 0))
    end
    
    # constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
        push!(pg_cons, @constraint(model, -kcl[gen["gen_bus"]] + pgmin[i] - pgmax[i] == 0))
    end

    # constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        push!(p_cons, @constraint(model, kcl[i] - kcl[j] + tmin[l] - tmax[l] + dclb[l] - dcub[l] == 0))
    end

    # constraints corresponding to primal variable ld
    for (i, bus) in ref[:bus]
        push!(ld_cons, @constraint(model, -bus["pd"]*kcl[i] - loadshed[i] <= bus["pd"]))
    end


    # primal objective == dual 
    @expression(model, primalobj_expr, sum( ref[:bus][i]["pd"] * ld[i] for i in keys(ref[:bus]) ) )
    @expression(model, dualobj_expr, 
                                sum( -ref[:bus][i]["pd"] * kcl[i] for i in keys(ref[:bus]) ) +
                                sum( ref[:gen][i]["pmin"] * pgmin[i] - ref[:gen][i]["pmax"] * pgmax[i] for i in keys(ref[:gen]) ) +
                                sum( -ref[:branch][l]["rate_a"] * tmin[l] - ref[:branch][l]["rate_a"] * tmax[l] for l in keys(ref[:branch]) ) +
                                sum( 0 + 0 for i in keys(ref[:branch]) ) + 
                                sum( -loadshed[i] for i in keys(ref[:bus]) )
                                )
    
    push!(pdeq_cons, @constraint(model, primalobj_expr - dualobj_expr == 0))





    return model

end
