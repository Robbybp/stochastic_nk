
function Lshaped(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, A, sense, l, u, master, solver)
    
    numscenarios = config["batchsize"]
    
    θ = getindex(master, :θ)
    x = getindex(master, :x)
    y = getindex(master, :y)
    numsubproblem_constr = length(sense)

    @objective(master, Max, sum(θ))
    

    function bendercuts(cb)
        xval = getvalue(x)
        yval = getvalue(y)
        θval = getvalue(θ)

        master_obj = sum(θval)

        b, c = create_vectors(scenarios, ref, config, Model(), xval, yval, numsubproblem_constr)

        @assert length(b) == numscenarios
        @assert length(c) == numscenarios

        subproblem_sol = Dict{Int,Any}()
        
        if config["parallel"] == "n"
            for s in 1:numscenarios
                subproblem_sol[s] = linprog(c[s], A, sense, b[s], l, u, solver)
                @assert subproblem_sol[s].status == :Optimal
            end
        else
            subproblem_sol = pmap((a1,a2,a3,a4,a5,a6,a7)->linprog(a1,a2,a3,a4,a5,a6,a7), 
                                  [c[s] for s in 1:numscenarios], 
                                  [A for s in 1:numscenarios], 
                                  [sense for s in 1:numscenarios], 
                                  [b[s] for s in 1:numscenarios], 
                                  [l for s in 1:numscenarios], 
                                  [u for s in 1:numscenarios], 
                                  [solver for s in 1:numscenarios])
        end

        
        # conversion of max problem to - min (-z) to use linprog
        subproblem_obj = [-subproblem_sol[s].objval for s in 1:numscenarios]

        if !isapprox(sum(subproblem_obj), master_obj, atol=1e-6)
            for s in 1:numscenarios
                @lazyconstraint(cb, θ[s] <= -sum( [master.ext[:proj_expr][s][i]*subproblem_sol[s].attrs[:lambda][i] for i in 1:numsubproblem_constr] ) )
            end
        end 

    end 
        
    addlazycallback(master, bendercuts)
    status = solve(master, suppress_warnings=true)
    
    obj = getobjectivevalue(master)/numscenarios

    sol = Dict{Symbol,Any}()
    sol[:x] = Int[]
    sol[:y] = Int[]
    
    for i in keys(ref[:nw][0][:branch])
        if getvalue(x[i]) > 0.9
            push!(sol[:x], i)
        end
    end

    for i in keys(ref[:nw][0][:gen])
        if getvalue(y[i]) > 0.9
            push!(sol[:y], i)
        end
    end

    config["final_mipgap"] = getobjgap(master)

    return status, obj, sol
end

function Lshaped_traditional(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, A, sense, l, u, master, solver)
    
    numscenarios = config["batchsize"]
    
    θ = getindex(master, :θ)
    x = getindex(master, :x)
    y = getindex(master, :y)
    numsubproblem_constr = length(sense)

    @objective(master, Max, sum(θ))
    
    lb = -1e5
    iteration = 1
    solve(master)
    ub = getobjectivevalue(master)

    while abs(ub-lb) > 1e-5*lb 
        xval = getvalue(x)
        yval = getvalue(y)

        b, c = create_vectors(scenarios, ref, config, Model(), xval, yval, numsubproblem_constr)

        @assert length(b) == numscenarios
        @assert length(c) == numscenarios

        subproblem_sol = Dict{Int,Any}()

        if config["parallel"] == "n"
            for s in 1:numscenarios
                subproblem_sol[s] = linprog(c[s], A, sense, b[s], l, u, solver)
                @assert subproblem_sol[s].status == :Optimal
            end
        else
            subproblem_sol = pmap((a1,a2,a3,a4,a5,a6,a7)->linprog(a1,a2,a3,a4,a5,a6,a7), 
                                  [c[s] for s in 1:numscenarios], 
                                  [A for s in 1:numscenarios], 
                                  [sense for s in 1:numscenarios], 
                                  [b[s] for s in 1:numscenarios], 
                                  [l for s in 1:numscenarios], 
                                  [u for s in 1:numscenarios], 
                                  [solver for s in 1:numscenarios])
        end
        
        # conversion of max problem to - min (-z) to use linprog
        subproblem_obj = [-subproblem_sol[s].objval for s in 1:numscenarios]
        lb = max(lb, sum(subproblem_obj))
        
        @constraint(master, [s=1:numscenarios], θ[s] <= -dot(master.ext[:proj_expr][s], subproblem_sol[s].attrs[:lambda]))
        solve(master)
        ub = getobjectivevalue(master)
        iteration += 1

        println(">> iteration: $iteration lb: $lb ub: $ub relgap: $(abs(ub-lb)/lb*100)")
    end 

           
    obj = getobjectivevalue(master)/numscenarios

    sol = Dict{Symbol,Any}()
    sol[:x] = Int[]
    sol[:y] = Int[]
    
    for i in keys(ref[:nw][0][:branch])
        if getvalue(x[i]) > 0.9
            push!(sol[:x], i)
        end
    end

    for i in keys(ref[:nw][0][:gen])
        if getvalue(y[i]) > 0.9
            push!(sol[:y], i)
        end
    end

    return :Optimal, obj, sol
end
