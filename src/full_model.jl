
function create_full_model(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}; model=Model())

    ref = ref[:nw][0]
    println(scenarios)
    numscenarios = config["batchsize"]
    config["budget"] == sum(scenarios[1,:])

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
    for k in 1:length(keys(ref[:gen]))
        gen_id = collect(keys(ref[:gen]))[k]
        gen = ref[:gen][gen_id]
        col_index = k
        @constraint(model, [s=1:numscenarios], pg[gen_id,s] >= (1-y[gen_id]*scenarios[s,col_index]) * gen["pmin"])
        @constraint(model, [s=1:numscenarios], pg[gen_id,s] <= (1-y[gen_id]*scenarios[s,col_index]) * gen["pmax"])
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

            #@constraint(model, p_fr + b*(va_fr - va_to) >= ref[:off_angmin]*x[i]*scenarios[s,col_index])
            #@constraint(model, p_fr + b*(va_fr - va_to) <= ref[:off_angmax]*x[i]*scenarios[s,col_index])
            @constraint(model, p_fr + b*(va_fr - va_to) >= -1000*x[i]*scenarios[s,col_index])
            @constraint(model, p_fr + b*(va_fr - va_to) <= 1000*x[i]*scenarios[s,col_index])
        end
    end

    # load shedding limit
    for (i, bus) in ref[:bus]
        @constraint(model, [s=1:numscenarios], ld[i,s] <= 1)
    end
    
    # dual constraints
    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        @constraint(model, [s=1:numscenarios], sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dclb[l,s] - dcub[l,s]) for (l,f,t) in bus_arcs) == 0)
    end

    # constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
         @constraint(model, [s=1:numscenarios], -kcl[gen["gen_bus"],s] + pgmin[i,s] - pgmax[i,s] == 0)
    end

    # constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        @constraint(model, [s=1:numscenarios], kcl[i,s] - kcl[j,s] + tmin[l,s] - tmax[l,s] + dclb[l,s] - dcub[l,s] == 0)
    end

    # constraints corresponding to primal variable ld
    for (i, bus) in ref[:bus]
        @constraint(model, [s=1:numscenarios], -bus["pd"]*kcl[i,s] - loadshed[i,s] <= bus["pd"])
    end

    # reformulation variables
    @variable(model, ypgmin[i in keys(ref[:gen]), s=1:numscenarios] >= 0)
    @variable(model, ypgmax[i in keys(ref[:gen]), s=1:numscenarios] >= 0)
    @variable(model, xtmin[i in keys(ref[:branch]), s=1:numscenarios] >= 0)
    @variable(model, xtmax[i in keys(ref[:branch]), s=1:numscenarios] >= 0)
    @variable(model, xdclb[i in keys(ref[:branch]), s=1:numscenarios] >= 0)
    @variable(model, xdcub[i in keys(ref[:branch]), s=1:numscenarios] >= 0)

    # reformulation constraints
    M = 1000
    for (i, gen) in ref[:gen]
        for s in 1:numscenarios
            add_reformulation(model, ypgmin[i,s], y[i], pgmin[i,s], M)
            add_reformulation(model, ypgmax[i,s], y[i], pgmax[i,s], M)
        end
    end

    for (l, branch) in ref[:branch]
        for s in 1:numscenarios
            add_reformulation(model, xtmin[l,s], x[l], tmin[l,s], M)
            add_reformulation(model, xtmax[l,s], x[l], tmax[l,s], M)
            add_reformulation(model, xdclb[l,s], x[l], dclb[l,s], M)
            add_reformulation(model, xdcub[l,s], x[l], dcub[l,s], M)
        end
    end
    

    # primal objective
    @expression(model, primalobj_expr[s=1:numscenarios], sum( ref[:bus][i]["pd"] * ld[i,s] for i in keys(ref[:bus]) ) )
    
    # dual objective
    kcl_expr = Any[]
    pgmin_expr = Any[]
    pgmax_expr = Any[]
    tmin_expr = Any[]
    tmax_expr = Any[]
    dc_expr = Any[]
    loadshed_expr = Any[]

    @expression(model, kcl_expr[s=1:numscenarios], sum( -ref[:bus][i]["pd"] * kcl[i,s] for i in keys(ref[:bus]) ) )
    
    @expression(model, loadshed_expr[s=1:numscenarios], sum( -loadshed[i,s] for i in keys(ref[:bus]) ) )

    @expression(model, pgmin_expr[s=1:numscenarios], 
                sum( ref[:gen][i]["pmin"] * pgmin[i,s] for i in keys(ref[:gen]) ) -
                sum( scenarios[s,i] * ref[:gen][collect(keys(ref[:gen]))[i]]["pmin"] * ypgmin[collect(keys(ref[:gen]))[i],s] for i in 1:length(keys(ref[:gen])) )
                )

    @expression(model, pgmax_expr[s=1:numscenarios], 
                sum( -ref[:gen][i]["pmax"] * pgmax[i,s] for i in keys(ref[:gen]) ) +
                sum( scenarios[s,i] * ref[:gen][collect(keys(ref[:gen]))[i]]["pmax"] * ypgmax[collect(keys(ref[:gen]))[i],s] for i in 1:length(keys(ref[:gen])) )
                )

    @expression(model, tmin_expr[s=1:numscenarios], 
                sum( -ref[:branch][l]["rate_a"] * tmin[l,s] for l in keys(ref[:branch]) ) +
                sum( ref[:branch][ref[:arcs_from][k][1]]["rate_a"] * scenarios[s,length(keys(ref[:gen]))+k] * xtmin[ref[:arcs_from][k][1],s] 
                    for k in 1:length(ref[:arcs_from]) )
                )

    @expression(model, tmax_expr[s=1:numscenarios], 
                sum( -ref[:branch][l]["rate_a"] * tmax[l,s] for l in keys(ref[:branch]) ) +
                sum( ref[:branch][ref[:arcs_from][k][1]]["rate_a"] * scenarios[s,length(keys(ref[:gen]))+k] * xtmax[ref[:arcs_from][k][1],s] 
                    for k in 1:length(ref[:arcs_from]) )
                )

    @expression(model, dc_expr[s=1:numscenarios], 
                sum( -1000 * scenarios[s,length(keys(ref[:gen]))+k] * xdclb[collect(keys(ref[:branch]))[k],s] for k in 1:length(ref[:branch]) ) +
                sum( -1000 * scenarios[s,length(keys(ref[:gen]))+k] * xdcub[collect(keys(ref[:branch]))[k],s] for k in 1:length(ref[:branch]) )
               )

    
    # primal objective == dual objective
    @constraint(model, [s=1:numscenarios], primalobj_expr[s] == kcl_expr[s] + pgmin_expr[s] + pgmax_expr[s] + tmin_expr[s] + tmax_expr[s] + dc_expr[s] + loadshed_expr[s])
    
    @objective(model, Max, sum(primalobj_expr)/numscenarios)
    return model

end

function add_reformulation(model, xy, x_bin, y_cont, y_ub)
    
    @constraint(model, xy >= y_cont - (1-x_bin)*y_ub)
    @constraint(model, xy <= y_cont)
    @constraint(model, xy <= y_ub*x_bin)

    return

end
