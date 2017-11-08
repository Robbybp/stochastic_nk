@everywhere using MathProgBase

function Lshaped_lazy(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, A, sense, l, u, master, solver)
    
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

        for s in 1:numscenarios
            subproblem_sol[s] = linprog(c[s], A, sense, b[s], l, u, solver)
            @assert subproblem_sol[s].status == :Optimal
        end
        
        # conversion of max problem to - min (-z) to use linprog
        subproblem_obj = [-subproblem_sol[s].objval for s in 1:numscenarios]

        if !isapprox(sum(subproblem_obj), master_obj, atol=1e-6)
            for s in 1:numscenarios
                @lazyconstraint(cb, θ[s] <= -dot(master.ext[:proj_expr][s], subproblem_sol[s].attrs[:lambda]))
            end
        end 

    end 
        
    addlazycallback(master, bendercuts)
    status = solve(master)
    
    obj = getobjectivevalue(master)/numscenarios

    sol = Dict{Symbol,Any}()
    sol[:x] = getvalue(x)
    sol[:y] = getvalue(y)

    return status, obj, sol
end

function Lshaped(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, A, sense, l, u, master, solver)
    
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

        for s in 1:numscenarios
            subproblem_sol[s] = linprog(c[s], A, sense, b[s], l, u, solver)
            @assert subproblem_sol[s].status == :Optimal
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
    sol[:x] = getvalue(x)
    sol[:y] = getvalue(y)

    return :Optimal, obj, sol
end
