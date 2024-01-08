""" run algorithm for stochastic N-k """
function run_stochastic(config::Dict, mp_file::String, scenario_file::String)
    data = PowerModels.parse_file(mp_file; validate=false)
    PowerModels.make_per_unit!(data)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    scenario_data = JSON.parsefile(scenario_file)
    num_scenarios = length(scenario_data)
    if (config["maximum_scenarios"] > num_scenarios)
        @warn "maximum scenarios requested $(config["maximum_scenarios"]) > number of scenarios in scenario file ($num_scenarios)"
        @warn "using $num_scenarios instead"
    end 
    filter!(p -> parse(Int64, first(p)) <= config["maximum_scenarios"], scenario_data)

    return solve_stochastic(config, data, ref, scenario_data)
end 

""" solve with lazy constraint callback """
function solve_stochastic(config::Dict, data::Dict, ref::Dict, scenarios::Dict)
    num_scenarios = length(scenarios)
    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # set_attribute(model, "LogToConsole", 0)
    set_attribute(model, "TimeLimit", config["timeout"])
    MOI.set(model, MOI.RelativeGapTolerance(), config["optimality_gap"] / 100.0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta[j in keys(scenarios)] <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])
    if config["use_separate_budgets"]
        @constraint(model, sum(x_line) == config["line_budget"])
        @constraint(model, sum(x_gen) == config["generator_budget"])
    end 

    # objective (expected load shed's SAA approximation)
    @objective(model, Max, sum(eta)/num_scenarios)

    TOL = 1E-6

    function inner_problem(cb_data)
        status = callback_node_status(cb_data, model)
        (status != MOI.CALLBACK_NODE_STATUS_INTEGER) && (return)
        current_x_line = Dict(i => JuMP.callback_value(cb_data, x_line[i]) for i in keys(ref[:branch]))
        current_x_gen = Dict(i => JuMP.callback_value(cb_data, x_gen[i]) for i in keys(ref[:gen]))
        current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
        current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect

        lck = Threads.ReentrantLock();
        Threads.@threads for val in collect(scenarios)
            s = first(val)
            scenario = last(val)
            cut_info = get_inner_solution(data, ref, current_gens, current_lines, scenario["gen"], scenario["branch"])
            woods_cut = @build_constraint(eta[s] <= cut_info.load_shed + sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]) + 
                sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)]))
            Threads.lock(lck) do 
                MOI.submit(model, MOI.LazyConstraint(cb_data), woods_cut) 
            end
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