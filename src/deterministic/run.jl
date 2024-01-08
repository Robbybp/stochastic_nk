"""
 * Implements the cutting plane algorithm in the paper
 *
 * Worst-Case Interdiction Analysis of Large-Scale Electric Power Grids 
 * DOI: 10.1109/TPWRS.2008.2004825
"""

""" run algorithm for determinsitic N-k """
function run_deterministic(config::Dict, mp_file::String)
    data = PowerModels.parse_file(mp_file; validate=false)
    PowerModels.make_per_unit!(data)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    return solve_deterministic(config, data, ref)
end 

""" solve with lazy constraint callback """
function solve_deterministic(config::Dict, data::Dict, ref::Dict)
    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # set_attribute(model, "LogToConsole", 0)
    set_attribute(model, "TimeLimit", config["timeout"])
    MOI.set(model, MOI.RelativeGapTolerance(), config["optimality_gap"] / 100.0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])
    if config["use_separate_budgets"]
        @constraint(model, sum(x_line) == config["line_budget"])
        @constraint(model, sum(x_gen) == config["generator_budget"])
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
        cut_info = get_inner_solution(data, ref, current_gens, current_lines)
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
    incumbent = SolutionDeterministic(current_lines, current_gens, objective_value, Dict())

    return ResultsDeterministic(
        iterations, objective_value, bound, run_time, rel_gap, incumbent
    )
end 