@everywhere using MathProgBase

function Lshaped_pmap(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, A, sense, l, u, master, solver, xval, yval)
    
    numscenarios = config["batchsize"]
    iteration_count = 1
    lb = -1e7
    ub = 1e7
    
    while abs(ub-lb)/abs(lb) > 1e-5
        
        θ = getindex(master, :θ)
        x = getindex(master, :x)
        y = getindex(master, :y)
        numsubproblem_constr = length(sense)
        
        b, c = create_vectors(scenarios, ref, config, Model(), xval, yval, numsubproblem_constr)

        @assert length(b) == numscenarios
        @assert length(c) == numscenarios
        
        subproblem_sol = Dict{Int,Any}()
        
        for s in 1:numscenarios
            subproblem_sol[s] = linprog(c[s], A, sense, b[s], l, u, solver)
            @assert subproblem_sol[s].status == :Optimal
        end

        subproblem_obj = [-subproblem_sol[s].objval for s in 1:numscenarios]
        if iteration_count > 1
            lb = max(sum(subproblem_obj)/numscenarios, lb)
        end
        
        @constraint(master, [s=1:numscenarios], θ[s] <= subproblem_obj[s] - dot(master.ext[:proj_expr][s], subproblem_sol[s].attrs[:lambda]))
        
        x_index = []
        y_index = []
        
        x_val = i -> getvalue(x[i])
        y_val = i -> getvalue(y[i])

        current_xvals = Dict( i => x_val(i) for i in keys(ref[:nw][0][:branch]) )
        x_index = collect(keys(filter( (i, val) -> val > 0.9, current_xvals)))

        current_yvals = Dict( i => y_val(i) for i in keys(ref[:nw][0][:gen]) )
        y_index = collect(keys(filter( (i, val) -> val > 0.9, current_xvals)))

        @constraint(master, sum(x[i] for i in x_index) + sum(y[i] for i in y_index) <= config["budget"] - 1)

        @objective(master, Max, sum(θ)/numscenarios)

        solve(master)
        ub = getobjectivevalue(master)

        println(">> iteration $iteration_count upperbound: $ub lowerbound: $lb")

        xval = getvalue(x)
        yval = getvalue(y)

        iteration_count += 1
        
    end

    obj = lowerbound
    sol[:x] = getvalue(getindex(master, :x))
    sol[:y] = getvalue(getindex(master, :y))

    return :Optimal, obj, sol
end
