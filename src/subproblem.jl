
function create_matrices(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, model=Model())

    ref = ref[:nw][0]
    numscenarios = config["batchsize"]
    M = 1000

    xvals = Dict([ i => 1 for i in keys(ref[:branch]) ])
    yvals = Dict([ i => 1 for i in keys(ref[:gen]) ])
    
    # dual variables
    @variable(model, kcl[i in keys(ref[:bus])])
    @variable(model, pgmin[i in keys(ref[:gen])]>= 0)
    @variable(model, pgmax[i in keys(ref[:gen])] >= 0)
    @variable(model, tmin[i in keys(ref[:branch])] >= 0)
    @variable(model, tmax[i in keys(ref[:branch])] >= 0)
    @variable(model, dclb[i in keys(ref[:branch])] >= 0)
    @variable(model, dcub[i in keys(ref[:branch])] >= 0)
    @variable(model, vamin[i in keys(ref[:branch])] >= 0)
    @variable(model, vamax[i in keys(ref[:branch])] >= 0)
    @variable(model, loadshed[i in keys(ref[:bus])] >= 0)

    # linearization variables
    @variable(model, ypgmax[i in keys(ref[:gen])] >= 0)
    @variable(model, xtmin[i in keys(ref[:branch])] >= 0)
    @variable(model, xtmax[i in keys(ref[:branch])] >= 0)
    @variable(model, xdclb[i in keys(ref[:branch])] >= 0)
    @variable(model, xdcub[i in keys(ref[:branch])] >= 0)
    @variable(model, xvamin[i in keys(ref[:branch])] >= 0)
    @variable(model, xvamax[i in keys(ref[:branch])] >= 0)

    # auxiliary definitions for formulating the constraints
    p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
    p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))

    # (a) constraints corresponding to primal variable va
    for (i, bus) in ref[:bus]
        bus_arcs = ref[:bus_arcs][i]
        @constraint(model, sum(PMs.calc_branch_y(ref[:branch][l])[2] * p_mag[(l,f,t)] * (dclb[l] - dcub[l]) + p_mag[(l,f,t)] * (vamin[l] - vamax[l]) for (l,f,t) in bus_arcs) == 0)
    end

    # (b) constraints corresponding to primal variable pg
    for (i, gen) in ref[:gen]
         @constraint(model, -kcl[gen["gen_bus"]] + pgmin[i] - pgmax[i] == 0)
    end

    # (c) constraints corresponding to primal variable p
    for (l,i,j) in ref[:arcs_from]
        @constraint(model, kcl[i] - kcl[j] + tmin[l] - tmax[l] + dclb[l] - dcub[l] == 0)
    end

    # (d) constraints corresponding to primal variable ld
    for (i, bus) in ref[:bus]
        @constraint(model, -bus["pd"]*kcl[i] - loadshed[i] <= bus["pd"])
    end

    # reformulation constraints
    for (i, gen) in ref[:gen]
        add_reformulation(model, ypgmax[i], yvals[i], pgmax[i], M)
    end

    for (l, branch) in ref[:branch]
        add_reformulation(model, xtmin[l], xvals[l], tmin[l], M)
        add_reformulation(model, xtmax[l], xvals[l], tmax[l], M)
        add_reformulation(model, xdclb[l], xvals[l], dclb[l], M)
        add_reformulation(model, xdcub[l], xvals[l], dcub[l], M)
        add_reformulation(model, xvamin[l], xvals[l], vamin[l], M)
        add_reformulation(model, xvamax[l], xvals[l], vamax[l], M)
    end
    
    sense = Char[]
    for i in 1:length(model.linconstr)
        constr_sense = JuMP.sense(model.linconstr[i])
        if constr_sense == :(<=)
            push!(sense, '<')
        elseif constr_sense == :(>=)
            push!(sense, '>')
        else
            push!(sense, '=')
        end
    end

    return JuMP.prepConstrMatrix(model), sense, model.colLower, model.colUpper

end 

function create_vectors(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, model=Model(), xval=xval, yval=yval, b_len=b_len)

    ref = ref[:nw][0]
    numscenarios = config["batchsize"]
    b = Dict{Int,Vector{Float64}}()
    c = Dict{Int,Vector{Float64}}()
    
    for s in 1:numscenarios
        b[s] = Float64[]
        
        model = Model()

        # dual variables
        @variable(model, kcl[i in keys(ref[:bus])])
        @variable(model, pgmin[i in keys(ref[:gen])]>= 0)
        @variable(model, pgmax[i in keys(ref[:gen])] >= 0)
        @variable(model, tmin[i in keys(ref[:branch])] >= 0)
        @variable(model, tmax[i in keys(ref[:branch])] >= 0)
        @variable(model, dclb[i in keys(ref[:branch])] >= 0)
        @variable(model, dcub[i in keys(ref[:branch])] >= 0)
        @variable(model, vamin[i in keys(ref[:branch])] >= 0)
        @variable(model, vamax[i in keys(ref[:branch])] >= 0)
        @variable(model, loadshed[i in keys(ref[:bus])] >= 0)

        # linearization variables
        @variable(model, ypgmax[i in keys(ref[:gen])] >= 0)
        @variable(model, xtmin[i in keys(ref[:branch])] >= 0)
        @variable(model, xtmax[i in keys(ref[:branch])] >= 0)
        @variable(model, xdclb[i in keys(ref[:branch])] >= 0)
        @variable(model, xdcub[i in keys(ref[:branch])] >= 0)
        @variable(model, xvamin[i in keys(ref[:branch])] >= 0)
        @variable(model, xvamax[i in keys(ref[:branch])] >= 0)

        # auxiliary definitions for formulating the constraints
        p_mag = Dict([((l,i,j), 1.0) for (l,i,j) in ref[:arcs_from]])
        p_mag = merge(p_mag, Dict([((l,j,i), -1.0) for (l,i,j) in ref[:arcs_from]]))
    

        # (a) constraints corresponding to primal variable va
        append!(b[s], zeros(length(ref[:bus])))

        # (b) constraints corresponding to primal variable pg
        append!(b[s], zeros(length(ref[:gen])))

        # (c) constraints corresponding to primal variable p
        append!(b[s], zeros(length(ref[:arcs_from])))

        # (d) constraints corresponding to primal variable ld
        for (i, bus) in ref[:bus]
            push!(b[s], bus["pd"])
        end

        # reformulation constraints
        for (i, gen) in ref[:gen]
            append_rhs(b[s], ypgmax[i], yval[i], pgmax[i], config["bounds"][:dual_pgmax][s][i])
        end

        for (l, branch) in ref[:branch]
            append_rhs(b[s], xtmin[l], xval[l], tmin[l], config["bounds"][:dual_tmin][s][l])
            append_rhs(b[s], xtmax[l], xval[l], tmax[l], config["bounds"][:dual_tmax][s][l])
            append_rhs(b[s], xdclb[l], xval[l], dclb[l], config["bounds"][:dual_dclb][s][l])
            append_rhs(b[s], xdcub[l], xval[l], dcub[l], config["bounds"][:dual_dcub][s][l])
            append_rhs(b[s], xvamin[l], xval[l], vamin[l], config["bounds"][:dual_vamin][s][l])
            append_rhs(b[s], xvamax[l], xval[l], vamax[l], config["bounds"][:dual_vamax][s][l])
        end

        # dual objective
        kcl_expr = Any[]
        pgmax_expr = Any[]
        tmin_expr = Any[]
        tmax_expr = Any[]
        dc_expr = Any[]
        va_expr = Any[]
        loadshed_expr = Any[]

        @expression(model, kcl_expr, sum( -ref[:bus][i]["pd"] * kcl[i] for i in keys(ref[:bus]) ) )

        @expression(model, loadshed_expr, sum( -loadshed[i] for i in keys(ref[:bus]) ) )

        @expression(model, pgmax_expr, 
                    sum( -ref[:gen][i]["pmax"] * pgmax[i] for i in keys(ref[:gen]) ) +
                    sum( scenarios[s,i] * ref[:gen][collect(keys(ref[:gen]))[i]]["pmax"] * ypgmax[collect(keys(ref[:gen]))[i]] for i in 1:length(keys(ref[:gen])) )
                   )

        @expression(model, tmin_expr, 
                    sum( -ref[:branch][l]["rate_a"] * tmin[l] for l in keys(ref[:branch]) ) +
                    sum( ref[:branch][ref[:arcs_from][k][1]]["rate_a"] * scenarios[s,length(keys(ref[:gen]))+k] * xtmin[ref[:arcs_from][k][1]] 
                        for k in 1:length(ref[:arcs_from]) )
                   )

        @expression(model, tmax_expr, 
                    sum( -ref[:branch][l]["rate_a"] * tmax[l] for l in keys(ref[:branch]) ) +
                    sum( ref[:branch][ref[:arcs_from][k][1]]["rate_a"] * scenarios[s,length(keys(ref[:gen]))+k] * xtmax[ref[:arcs_from][k][1]] 
                        for k in 1:length(ref[:arcs_from]) )
                   )

        t_min = ref[:off_angmin]
        t_max = ref[:off_angmax]

        @expression(model, dc_expr, 
                    sum( - PMs.calc_branch_y(ref[:branch][collect(keys(ref[:branch]))[k]])[2] * t_min * scenarios[s,length(keys(ref[:gen]))+k] * xdclb[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) ) +
                    sum( PMs.calc_branch_y(ref[:branch][collect(keys(ref[:branch]))[k]])[2] * t_max * scenarios[s,length(keys(ref[:gen]))+k] * xdcub[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) )
                   )

        @expression(model, va_expr, 
                    sum( ref[:branch][collect(keys(ref[:branch]))[k]]["angmin"] * vamin[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) ) +
                    sum( (-ref[:branch][collect(keys(ref[:branch]))[k]]["angmin"] + t_min) * scenarios[s,length(keys(ref[:gen]))+k] * xvamin[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) ) -
                    sum( ref[:branch][collect(keys(ref[:branch]))[k]]["angmax"] * vamax[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) ) -
                    sum( (-ref[:branch][collect(keys(ref[:branch]))[k]]["angmax"] + t_max) * scenarios[s,length(keys(ref[:gen]))+k] * xvamax[collect(keys(ref[:branch]))[k]] for k in 1:length(ref[:branch]) )
                   )

        @expression(model, dualobj_expr, kcl_expr + pgmax_expr + tmin_expr + tmax_expr + dc_expr + va_expr + loadshed_expr)
        
        @objective(model, Min, -dualobj_expr)
        
        c[s] = JuMP.prepAffObjective(model)

        @assert length(b[s]) == b_len
        
        end

    return b, c
    
end


function add_reformulation(model, xy, x_bin, y_cont, y_ub)
    
    @constraint(model, xy >= y_cont - (1-x_bin)*y_ub)
    @constraint(model, xy <= y_cont)
    @constraint(model, xy <= y_ub*x_bin)

    return

end

function append_rhs(b, xy, x_bin, y_cont, y_ub)
    
    push!(b, (x_bin-1)*y_ub)
    push!(b, 0)
    push!(b, y_ub*x_bin)

    return

end
