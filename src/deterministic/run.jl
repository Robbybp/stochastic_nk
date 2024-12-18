"""
 * Implements the cutting plane algorithm in the paper
 *
 * Worst-Case Interdiction Analysis of Large-Scale Electric Power Grids 
 * DOI: 10.1109/TPWRS.2008.2004825
"""

""" run algorithm for determinsitic N-k """
function run_deterministic(cliargs::Dict, mp_file::String)::Results
    data = PowerModels.parse_file(mp_file; validate=false)
    PowerModels.make_per_unit!(data)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    return solve_deterministic(cliargs, data, ref)
end 

""" Method that accepts the PowerModels data dictionary directly """
function run_deterministic(cliargs::Dict, data::Dict{String,Any})::Results
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    return solve_deterministic(cliargs, data, ref)
end

""" solve with lazy constraint callback """
function solve_deterministic(cliargs::Dict, data::Dict, ref::Dict)::Results
    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    # set_attribute(model, "LogToConsole", 0)
    set_attribute(model, "TimeLimit", cliargs["timeout"])
    MOI.set(model, MOI.RelativeGapTolerance(), cliargs["optimality_gap"] / 100.0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    if cliargs["interdict_buses"]
        @variable(model, x_bus[i in keys(ref[:bus])], Bin)
    end

    # budget constraints 
    if cliargs["use_separate_budgets"]
        @constraint(model, sum(x_line) == cliargs["line_budget"])
        @constraint(model, sum(x_gen) == cliargs["generator_budget"])
    elseif cliargs["interdict_buses"]
        @constraint(model, sum(x_bus) == cliargs["budget"])
    else
        @constraint(model, sum(x_line) + sum(x_gen) == cliargs["budget"])
    end 

    # Logic constraints: If a bus is interdicted, incident generators and lines
    # are disrupted.
    if cliargs["interdict_buses"]
        @constraint(model,
            #[i in keys(ref[:bus]), (l, ibus, jbus) in ref[:bus_arcs][i]],
            [(l, ibus, jbus) in ref[:arcs]],
            x_line[l] <= x_bus[ibus] + x_bus[jbus]
        )
        @constraint(model,
            [(l, ibus, jbus) in ref[:arcs]],
            x_bus[ibus] <= x_line[l]
        )
        @constraint(model,
            [(l, ibus, jbus) in ref[:arcs]],
            x_bus[jbus] <= x_line[l]
        )
        @constraint(model,
            [i in keys(ref[:bus]), g in ref[:bus_gens][i]],
            x_bus[i] == x_gen[g]
        )
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
        if cliargs["interdict_buses"]
            # TODO: Could just define these variables and fix them to zero if we're
            # not interdicting buses. This should be profiled to make sure it's not a
            # performance hit.
            current_x_bus = Dict(i => JuMP.callback_value(cb_data, x_bus[i]) for i in keys(ref[:bus]))
            current_buses = filter!(z -> last(z) > TOL, current_x_bus) |> keys |> collect
        else
            current_buses = nothing
        end
        # TODO: Accept current bus status in get_inner_solution
        cut_info = get_inner_solution(
            data,
            ref,
            current_gens,
            current_lines;
            solver=cliargs["inner_solver"],
        )
        woods_cut = @build_constraint(
            eta <= (
                cut_info.load_shed
                + sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)])
                + sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)])
            )
        )
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

    if cliargs["interdict_buses"]
        current_x_bus = Dict(i => JuMP.value(x_bus[i]) for i in keys(ref[:bus]))
        current_buses = filter!(z -> last(z) > TOL, current_x_bus) |> keys |> collect
    else
        current_buses = []
    end

    # TODO: If we're interdicting buses, should we report only buses, or the lines and
    # generators as well? I think it makes sense to report all interdicted components.
    incumbent = Solution(current_lines, current_gens, current_buses, objective_value, Dict())

    return Results(
        iterations, objective_value, bound, run_time, rel_gap, incumbent
    )
end 
