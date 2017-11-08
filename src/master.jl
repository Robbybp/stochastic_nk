
function create_master_model(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, model=Model())
    
    ref = ref[:nw][0]
    numscenarios = config["batchsize"]

    # interdiction variables
    @variable(model, x[i in keys(ref[:branch])], Bin)
    @variable(model, y[i in keys(ref[:gen])], Bin)

    # lifted variables for multi-cut Lshaped
    @variable(model, θ[s=1:numscenarios] <= config["theta_ub"][s])

    # master constraints
    @constraint(model, sum(x) + sum(y) == config["budget"])
    (config["interdict"] == "l") && @constraint(model, sum(y) == 0)

    return model

end

function create_expression_vectors(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}; master_model=master)
    
    ref = ref[:nw][0]
    numscenarios = config["batchsize"]
    M = 1000

    # add references to tightened variable bounds to the configuration dictionary
    if !haskey(config, "bounds")
        config["bounds"] = Dict{Symbol,Any}()

        config["bounds"][:dual_pgmax] = Dict{Any,Any}()

        config["bounds"][:dual_tmin] = Dict{Any,Any}()
        config["bounds"][:dual_tmax] = Dict{Any,Any}()
        config["bounds"][:dual_dclb] = Dict{Any,Any}()
        config["bounds"][:dual_dcub] = Dict{Any,Any}()
        config["bounds"][:dual_vamin] = Dict{Any,Any}()
        config["bounds"][:dual_vamax] = Dict{Any,Any}()

        for s in 1:numscenarios
            config["bounds"][:dual_pgmax][s] = Dict{Any,Float64}( i => M for i in keys(ref[:gen]) )

            config["bounds"][:dual_tmin][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
            config["bounds"][:dual_tmax][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
            config["bounds"][:dual_dclb][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
            config["bounds"][:dual_dcub][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
            config["bounds"][:dual_vamin][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
            config["bounds"][:dual_vamax][s] = Dict{Any,Float64}( i => M for i in keys(ref[:branch]) )
        end
    end

    ext = Dict{Symbol,Any}()
    ext[:proj_expr] = Dict{Int,Any}()
    
    x = getindex(master_model, :x)
    y = getindex(master_model, :y)

    for s in 1:numscenarios
        ext[:proj_expr][s] = Any[]
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
        append!(ext[:proj_expr][s], zeros(length(ref[:bus])))

        # (b) constraints corresponding to primal variable pg
        append!(ext[:proj_expr][s], zeros(length(ref[:gen])))

        # (c) constraints corresponding to primal variable p
        append!(ext[:proj_expr][s], zeros(length(ref[:arcs_from])))

        # (d) constraints corresponding to primal variable ld
        for (i, bus) in ref[:bus]
            push!(ext[:proj_expr][s], bus["pd"])
        end

        # reformulation constraints
        for (i, gen) in ref[:gen]
            append_rhsx(ext[:proj_expr][s], ypgmax[i], y[i], pgmax[i], config["bounds"][:dual_pgmax][s][i])
        end

        for (l, branch) in ref[:branch]
            append_rhsx(ext[:proj_expr][s], xtmin[l], x[l], tmin[l], config["bounds"][:dual_tmin][s][l])
            append_rhsx(ext[:proj_expr][s], xtmax[l], x[l], tmax[l], config["bounds"][:dual_tmax][s][l])
            append_rhsx(ext[:proj_expr][s], xdclb[l], x[l], dclb[l], config["bounds"][:dual_dclb][s][l])
            append_rhsx(ext[:proj_expr][s], xdcub[l], x[l], dcub[l], config["bounds"][:dual_dcub][s][l])
            append_rhsx(ext[:proj_expr][s], xvamin[l], x[l], vamin[l], config["bounds"][:dual_vamin][s][l])
            append_rhsx(ext[:proj_expr][s], xvamax[l], x[l], vamax[l], config["bounds"][:dual_vamax][s][l])
        end        

    end
    
    return ext

end

function append_rhsx(bx, xy, x_bin, y_cont, y_ub)
    
    push!(bx, (x_bin-1)*y_ub)
    push!(bx, 0)
    push!(bx, y_ub*x_bin)

    return

end
