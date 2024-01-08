""" run algorithm for stochastic N-k """
function run_stochastic(config::Dict, mp_file::String, scenario_file::String)
    data = PowerModels.parse_file(mp_file)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    scendatatemp = JSON.parsefile(scenario_file)
    if config["maximum_scenarios"]<length(scendatatemp)
        scendatatemp = JSON.parsefile(scenario_file)
        scendata=Dict()
        counter=0
        for (s,scen) in scendatatemp
            if counter<=config["maximum_scenarios"]
                scendata[s] = scen
                counter+=1
            else
                break
            end
        end
    else
        scendata = scendatatemp
    end

    if (config["use_iterative"])
        fields, field_chars = get_table_config(problem = :deterministic)
        print_table_header(fields, field_chars)
        
        start_time = now() 

        return solve_stochastic_iterative(
            config, data, ref, scendata;
            start_time = start_time, 
            fields = fields, field_chars = field_chars)
    end  
    
    return solve_stochastic_lazy(config, data, ref, scendata)
end 


""" iterative solver without any lazy callbacks """
function solve_stochastic_iterative(config::Dict, data::Dict, ref::Dict, scendata::Dict; 
    start_time = now(), fields = nothing, field_chars = nothing)
    if isnothing(fields) || isnothing(field_chars)
        fields, field_chars = get_table_config(problem = :deterministic)
    end 

    lb = -Inf
    ub = Inf
    incumbent = SolutionStochastic()
    TOL = 1E-6
    rel_gap = Inf

    print_table_cp(0, lb, ub, rel_gap, start_time, fields, field_chars)

    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    set_attribute(model, "LogToConsole", 0)

    # eta represents the min load shed 
    @variable(model, 0 <= eta[j in keys(scendata)] <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) == config["line_budget"])
    @constraint(model, sum(x_gen) == config["generator_budget"])
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])

    # objective 
    @objective(model, Max, sum(eta)/length(scendata))

    get_ub(model) = JuMP.objective_value(model)

    iteration = 0

    cut_info = Dict()
    while (rel_gap > config["optimality_gap"])
        iteration += 1
        JuMP.optimize!(model)
        (get_ub(model) < ub) && (ub = get_ub(model))
        current_x_line = Dict(i => JuMP.value(x_line[i]) for i in keys(ref[:branch]))
        current_x_gen = Dict(i => JuMP.value(x_gen[i]) for i in keys(ref[:gen]))
        current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
        current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
        for (s, scen) in scendata
            current_lines_temp = copy(current_lines)
            for b in scen["branch"]
                if !(b in current_lines_temp)
                    append!(current_lines_temp,b)
                end
            end
            current_gens_temp = copy(current_gens)
            if "gen" in keys(scen)
                for g in scen["gen"]
                    if !(g in current_gens_temp)
                        append!(current_gens_temp,g)
                    end
                end
            end
            cut_info[s] = get_inner_solution(data, ref, current_gens_temp, current_lines_temp)
        end
        
        

        # update lower bound
        inner_objective = sum(cut_info[s].load_shed for s in keys(scendata))/length(scendata) 
        if inner_objective > lb 
            lb = inner_objective
            incumbent = SolutionStochastic(current_lines, current_gens, inner_objective, Dict())
        end 

        # add solution removal cut 
        if (config["generator_budget"] == 0)
            @constraint(model, sum([x_line[i] for i in current_lines]) <= config["budget"] - 1)
        else 
            @constraint(model, sum([x_line[i] for i in current_lines]) + sum([x_gen[i] for i in current_gens]) <= config["budget"] - 1)
        end 
    
        # add cutting planes   (inner objective specific to each scenario?)
        for (s, scen) in scendata
            @constraint(model, eta[s] <= cut_info[s].load_shed + 
                sum([cut_info[s].pg[i] * x_gen[i] for i in keys(cut_info[s].pg)]) + 
                sum([cut_info[s].p[i] * x_line[i] for i in keys(cut_info[s].p)])
            )
        end

        rel_gap = (ub - lb)/(1E-6 + abs(lb)) * 100
        (rel_gap < 0.0) && (ub = lb; rel_gap = 0.0)
        print_table_cp(iteration, lb, ub, rel_gap, start_time, fields, field_chars)
        
        if time_limit_reached(start_time, float(config["timeout"]))
            break
        end 
    end 

    print_table_footer(fields, field_chars)
    return ResultsStochastic(
        iteration, lb, ub, get_time(start_time), rel_gap, incumbent
    )
end 



""" solve with lazy constraint callback """
function solve_stochastic_lazy(config::Dict, data::Dict, ref::Dict, scendata::Dict)
    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # set_attribute(model, "LogToConsole", 0)
    set_attribute(model, "TimeLimit", config["timeout"])
    MOI.set(model, MOI.RelativeGapTolerance(), config["optimality_gap"] / 100.0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta[j in keys(scendata)] <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) == config["line_budget"])
    @constraint(model, sum(x_gen) == config["generator_budget"])
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])

    # objective 
    @objective(model, Max, sum(eta)/length(scendata))

    TOL = 1E-6

    function inner_problem(cb_data)
        status = callback_node_status(cb_data, model)
        (status != MOI.CALLBACK_NODE_STATUS_INTEGER) && (return)
        current_x_line = Dict(i => JuMP.callback_value(cb_data, x_line[i]) for i in keys(ref[:branch]))
        current_x_gen = Dict(i => JuMP.callback_value(cb_data, x_gen[i]) for i in keys(ref[:gen]))
        current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
        current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
        for (s, scen) in scendata
            current_lines_temp = copy(current_lines)
            for b in scen["branch"]
                if !(b in current_lines_temp)
                    append!(current_lines_temp,b)
                end
            end
            current_gens_temp = copy(current_gens)
            if "gen" in keys(scen)
                for g in scen["gen"]
                    if !(g in current_gens_temp)
                        append!(current_gens_temp,g)
                    end
                end
            end
            cut_info = get_inner_solution(data, ref, current_gens_temp, current_lines_temp)
            woods_cut = @build_constraint(eta[s] <= cut_info.load_shed + sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]) + 
                sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)]))
            MOI.submit(model, MOI.LazyConstraint(cb_data), woods_cut)
        end
    end

    MOI.set(model, MOI.LazyConstraintCallback(), inner_problem)
    JuMP.optimize!(model)

    iterations = 0
    run_time = JuMP.solve_time(model)
    objective_value = JuMP.objective_value(model)
    bound = JuMP.objective_bound(model)
    rel_gap = JuMP.relative_gap(model)
    current_x_line = Dict(i => JuMP.value(x_line[i]) for i in keys(ref[:branch]))
    current_x_gen = Dict(i => JuMP.value(x_gen[i]) for i in keys(ref[:gen]))
    current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
    current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
    incumbent = SolutionStochastic(current_lines, current_gens, objective_value, Dict())

    return ResultsStochastic(
        iterations, objective_value, bound, run_time, rel_gap, incumbent
    )
end 


# The functions "get_inner_solution" and "get_inner_solution" are defined in "deterministic/run.jl"