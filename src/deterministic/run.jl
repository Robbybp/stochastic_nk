"""
 * Implements the cutting plane algorithm in the paper
 *
 * Worst-Case Interdiction Analysis of Large-Scale Electric Power Grids 
 * DOI: 10.1109/TPWRS.2008.2004825
"""

const GRB_ENV = Gurobi.Env()

""" run algorithm for determinsitic N-k """
function run_deterministic(config::Dict, mp_file::String)
    data = PowerModels.parse_file(mp_file)
    add_total_load_info(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]

    if (config["use_lazy"] == false)
        fields, field_chars = get_table_config(problem = :deterministic)
        print_table_header(fields, field_chars)
        
        start_time = now() 

        return solve_deterministic_iterative(
            config, data, ref;
            start_time = start_time, 
            fields = fields, field_chars = field_chars)
    end  
    
    return solve_deterministic_lazy(config, data, ref)
end 

""" iterative solver without any lazy callbacks """
function solve_deterministic_iterative(config::Dict, data::Dict, ref::Dict; 
    start_time = now(), fields = nothing, field_chars = nothing)
    if isnothing(fields) || isnothing(field_chars)
        fields, field_chars = get_table_config(problem = :deterministic)
    end 

    lb = -Inf
    ub = Inf
    incumbent = SolutionDeterministic()
    TOL = 1E-6
    rel_gap = Inf

    print_table_cp(0, lb, ub, rel_gap, start_time, fields, field_chars)

    model = direct_model(Gurobi.Optimizer(GRB_ENV))
    set_attribute(model, "LogToConsole", 0)

    # eta represents the min load shed 
    @variable(model, 1E-6 <= eta <= 1E6)
    # interdiction variables
    @variable(model, x_line[i in keys(ref[:branch])], Bin)
    @variable(model, x_gen[i in keys(ref[:gen])], Bin)
    
    # budget constraints 
    @constraint(model, sum(x_line) == config["line_budget"])
    @constraint(model, sum(x_gen) == config["generator_budget"])
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])

    # objective 
    @objective(model, Max, eta)

    get_ub(model) = JuMP.objective_value(model)

    iteration = 0

    while (rel_gap > config["optimality_gap"])
        iteration += 1
        JuMP.optimize!(model)
        (get_ub(model) < ub) && (ub = get_ub(model))
        current_x_line = Dict(i => JuMP.value(x_line[i]) for i in keys(ref[:branch]))
        current_x_gen = Dict(i => JuMP.value(x_gen[i]) for i in keys(ref[:gen]))
        current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
        current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
        cut_info = get_inner_solution(data, ref, current_gens, current_lines)

        # update lower bound
        inner_objective = cut_info.load_shed 
        if inner_objective > lb 
            lb = inner_objective
            incumbent = SolutionDeterministic(current_lines, current_gens, inner_objective, Dict())
        end 

        # add solution removal cut 
        if (config["generator_budget"] == 0)
            @constraint(model, sum([x_line[i] for i in current_lines]) <= config["budget"] - 1)
        else 
            @constraint(model, sum([x_line[i] for i in current_lines]) + 
                sum([x_gen[i] for i in current_gens]) <= config["budget"] - 1)
        end 
    
        # add cutting plane 
        @constraint(model, eta <= inner_objective + 
            sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]) + 
            sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)])
        )

        rel_gap = (ub - lb)/(1E-6 + abs(lb)) * 100
        (rel_gap < 0.0) && (ub = lb; rel_gap = 0.0)
        print_table_cp(iteration, lb, ub, rel_gap, start_time, fields, field_chars)
        
        if time_limit_reached(start_time, float(config["timeout"]))
            break
        end 
    end 

    print_table_footer(fields, field_chars)
    return ResultsDeterministic(
        iteration, lb, ub, get_time(start_time), rel_gap, incumbent
    )
end 

""" solve with lazy constraint callback """
function solve_deterministic_lazy(config::Dict, data::Dict, ref::Dict)
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
    @constraint(model, sum(x_line) == config["line_budget"])
    @constraint(model, sum(x_gen) == config["generator_budget"])
    @constraint(model, sum(x_line) + sum(x_gen) == config["budget"])

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

""" get load shed and power flow solution on interdictable components""" 
function get_inner_solution(data, ref, generators::Vector, lines::Vector)::NamedTuple
    case_data = data
    case = deepcopy(case_data)
    for i in generators 
        case["gen"][string(i)]["gen_status"] = 0
    end 

    for i in lines 
        case["branch"][string(i)]["br_status"] = 0
    end 

    PowerModels.propagate_topology_status!(case)
    lp_optimizer = JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), "LogToConsole" => 0
    )
    
    pm = instantiate_model(case, DCPPowerModel, PowerModels._build_mld)
    result = optimize_model!(pm, optimizer = lp_optimizer)

    # result = solve_model(case, DCPPowerModel, lp_optimizer, PowerModels._build_mld)
    load_served = [load["pd"] for (_, load) in result["solution"]["load"]] |> sum
    load_shed = case_data["total_load"] - load_served

    pg = Dict(i => result["solution"]["gen"][string(i)]["pg"]
        for i in keys(ref[:gen]) if haskey(result["solution"]["gen"], string(i)))
    
    p = Dict(i => max(
        abs(result["solution"]["branch"][string(i)]["pf"]), 
        abs(result["solution"]["branch"][string(i)]["pt"])
        ) for i in keys(ref[:branch]) if haskey(result["solution"]["branch"], string(i)) 
    )

    return (load_shed = load_shed, pg = pg, p = p)
end 

""" get load shed and power flow solution on fractional interdictable components""" 
function get_inner_solution(data, ref, generators::Dict{Int,Float64}, lines::Dict{Int,Float64})::NamedTuple
    case_data = data
    case = deepcopy(case_data)
    for i in generators 
        case["gen"][string(i)]["gen_status"] = 0
    end 

    for i in lines 
        case["branch"][string(i)]["br_status"] = 0
    end 

    PowerModels.propagate_topology_status!(case)
    lp_optimizer = JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), "LogToConsole" => 0
    )
    
    pm = instantiate_model(case, DCPPowerModel, PowerModels._build_mld)
    result = optimize_model!(pm, optimizer = lp_optimizer)

    # result = solve_model(case, DCPPowerModel, lp_optimizer, PowerModels._build_mld)
    load_served = [load["pd"] for (_, load) in result["solution"]["load"]] |> sum
    load_shed = case_data["total_load"] - load_served

    pg = Dict(i => result["solution"]["gen"][string(i)]["pg"]
        for i in keys(ref[:gen]) if haskey(result["solution"]["gen"], string(i)))
    
    p = Dict(i => max(
        abs(result["solution"]["branch"][string(i)]["pf"]), 
        abs(result["solution"]["branch"][string(i)]["pt"])
        ) for i in keys(ref[:branch]) if haskey(result["solution"]["branch"], string(i)) 
    )

    return (load_shed = load_shed, pg = pg, p = p)
end 
