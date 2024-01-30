""" run algorithm for EVP N-k """
function run_evp(cliargs::Dict, mp_file::String, scenario_file::String)::Results
    data = PowerModels.parse_file(mp_file; validate=false)
    PowerModels.make_per_unit!(data)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    scenario_data = JSON.parsefile(scenario_file)
    num_scenarios = length(scenario_data)
    if (cliargs["maximum_scenarios"] > num_scenarios)
        @warn "maximum scenarios requested $(cliargs["maximum_scenarios"]) > number of scenarios in scenario file ($num_scenarios)"
        @warn "using $num_scenarios instead"
    end 
    filter!(p -> parse(Int64, first(p)) <= cliargs["maximum_scenarios"], scenario_data)

    scenario_generators = Dict{Int,Float64}()
    scenario_lines = Dict{Int,Float64}()
    for (_, scenario) in scenario_data
        lines = scenario["branch"]
        gens = scenario["gen"]
        for i in gens 
            (haskey(scenario_generators, i)) && (scenario_generators[i] += 1.0; continue)
            scenario_generators[i] = 1.0  
        end 
        for i in lines 
            (haskey(scenario_lines, i)) && (scenario_lines[i] += 1.0; continue)
            scenario_lines[i] = 1.0 
        end 
    end 
    for (i, val) in scenario_generators 
        scenario_generators[i] = val/length(scenario_data) 
    end 
    for (i, val) in scenario_lines 
        scenario_lines[i] = val/length(scenario_data) 
    end
    return solve_evp(cliargs, data, ref, scenario_generators, scenario_lines)
end 

""" solve with lazy constraint callback """
function solve_evp(cliargs::Dict, 
    data::Dict, ref::Dict, 
    scenario_generators::Dict,
    scenario_lines::Dict)::Results
    
    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # set_attribute(model, "LogToConsole", 0)
    set_attribute(model, "TimeLimit", cliargs["timeout"])
    MOI.set(model, MOI.RelativeGapTolerance(), cliargs["optimality_gap"] / 100.0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) + sum(x_gen) == cliargs["budget"])
    if cliargs["use_separate_budgets"]
        @constraint(model, sum(x_line) == cliargs["line_budget"])
        @constraint(model, sum(x_gen) == cliargs["generator_budget"])
    end 

    # objective
    @objective(model, Max, eta)

    TOL = 1E-6

    function inner_problem(cb_data)
        status = callback_node_status(cb_data, model)
        (status != MOI.CALLBACK_NODE_STATUS_INTEGER) && (return)
        current_x_line = Dict(i => JuMP.callback_value(cb_data, x_line[i]) for i in keys(ref[:branch]))
        current_x_gen = Dict(i => JuMP.callback_value(cb_data, x_gen[i]) for i in keys(ref[:gen]))
        current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
        current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
        cut_info = get_inner_solution(data, ref, 
            Dict{Int,Float64}(i => 1.0 for i in current_gens), 
            Dict{Int,Float64}(i => 1.0 for i in current_lines), 
            scenario_generators, scenario_lines; 
            solver=cliargs["inner_solver"])
        woods_cut = @build_constraint(eta <= cut_info.load_shed + sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]) + 
            sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)]))
        MOI.submit(model, MOI.LazyConstraint(cb_data), woods_cut)
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
    incumbent = Solution(current_lines, current_gens, objective_value, Dict())

    return Results(
        iterations, objective_value, bound, run_time, rel_gap, incumbent
    )
end 