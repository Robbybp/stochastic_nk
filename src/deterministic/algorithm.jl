"""
 * Implements the cutting plane algorithm in the paper
 *
 * "Worst-Case Interdiction Analysis of Large-Scale Electric Power Grids" 
 * doi:10.1109/TPWRS.2008.2004825
"""

function deterministic_cutting_plane()

    function cutting_plane(input_data, input_params, 
        solver_params::Dict; ignore_separate_budgets = true)::Tuple{Solution, PriorityQueue}
    
        fields, field_chars = get_table_config(algo = :cp)
        print_table_header(fields, field_chars)
        
        solutions = PriorityQueue{Solution, Float64}()
        start_time = now() 
    
        solutions = cp(input_data, input_params, solver_params, 
            solutions = solutions, 
            ignore_separate_budgets = ignore_separate_budgets,
            start_time = start_time, 
            fields = fields, 
            field_chars = field_chars)
    
        while (length(solutions) < solver_params["max_num_solutions"])
            solutions = cp(input_data, input_params, solver_params, 
                solutions = solutions, 
                ignore_separate_budgets = ignore_separate_budgets, 
                start_time = start_time,
                fields = fields, 
                field_chars = field_chars)
            if time_limit_reached(start_time, solver_params["time_limit_in_seconds"])
                print_table_footer(fields, field_chars)
                return (get_best(solutions), solutions)
            end 
        end 
    
        print_table_footer(fields, field_chars)
        return (get_best(solutions), solutions)
    end 
    
    function cp(input_data, input_params,
        solver_params::Dict; 
        solutions = PriorityQueue{Solution, Float64}(),
        ignore_separate_budgets = true,
        start_time = now(), 
        fields = nothing, 
        field_chars = nothing
    )::PriorityQueue
    
        if isnothing(fields) || isnothing(field_chars)
            fields, field_chars = get_table_config(algo = :cp)
        end 
        
        if (isempty(solutions)) 
            x = construct_initial_solution(input_data, input_params)
            add_to_solutions(solutions, x, max_num_solutions = solver_params["max_num_solutions"])
        end 
    
        interdictable_components = input_data.interdictable_components 
        probabilities = input_data.probabilities 
        bus_prob = probabilities.bus_prob 
        line_prob = probabilities.line_prob 
        gen_prob = probabilities.gen_prob
        
        lb = -Inf
        ub = Inf
        TOL = 1E-6
    
    
        print_table_cp(solutions |> length, 0, 
            get_best(solutions).cost, get_worst(solutions).cost, start_time, lb, ub, Inf, fields, field_chars)
    
        model = JuMP.Model(lp_optimizer)
    
        @variable(model, 1E-6 <= eta <= 1E6)
        @variable(model, x_bus[i in interdictable_components.buses], Bin)
        @variable(model, x_line[i in interdictable_components.lines], Bin)
        @variable(model, x_gen[i in interdictable_components.gens], Bin)
        @variable(model, p <= 0)
        @variable(model, y <= 1E6)
        
        if ignore_separate_budgets == false
            @constraint(model, sum(x_bus) == input_params["bus_budget"])
            @constraint(model, sum(x_line) == input_params["line_budget"])
            @constraint(model, sum(x_gen) == input_params["generator_budget"])
        else 
            @constraint(model, sum(x_bus) + sum(x_line) + sum(x_gen) == input_params["total_budget"])
        end 
    
        # logical constraints
        for i in interdictable_components.buses 
            incident_lines = input_data.bus_data.lines[i]
            bus_gens = get(input_data.bus_data.gens, i, [])
            for (l, _, _) in incident_lines 
                if (l in interdictable_components.lines)
                    @constraint(model, x_line[l] <= 1 - x_bus[i])
                end 
            end 
            for g in bus_gens 
                if (g in interdictable_components.gens)
                    @constraint(model, x_gen[g] <= 1 - x_bus[i])
                end 
            end 
        end 
    
        @constraint(model, p == 
            sum(x_bus[i] * log(bus_prob[i]) for i in interdictable_components.buses; init = 0.0) + 
            sum(x_line[i] * log(line_prob[i]) for i in interdictable_components.lines; init = 0.0) + 
            sum(x_gen[i] * log(gen_prob[i]) for i in interdictable_components.gens; init = 0.0)
        )   
    
        for solution in solutions
            z = first(solution)
            @constraint(model, sum([x_bus[i] for i in z.buses]; init = 0.0) + 
                sum([x_line[i] for i in z.lines]; init = 0.0) + 
                sum([x_gen[i] for i in z.generators]; init = 0.0) <= input_params["total_budget"] - 1)
    
            # cut_info = get_solution_values(input_data, z.generators, z.buses, z.lines)
            # @constraint(model, eta <= cut_info.load_shed + 
            #     sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]; init = 0.0) + 
            #     sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)]; init = 0.0) +
            #     sum([cut_info.pb[i] * x_bus[i] for i in keys(cut_info.pb)]; init = 0.0)
            # )
        end 
    
        @objective(model, Max, y)
        # add outer approximation of constraint y ⩽ p + log η using lazy constraint callback:  y ⩽ p + log η₀ + 1/η₀(η - η₀)
        # outer approximation of y ⩽ f(x) at x = xᵏ : y ⩽ f(xᵏ) + ∇f(xᵏ)⋅(x-xᵏ)    
    
        function outer_approximate(cb_data)
            status = callback_node_status(cb_data, model)
            (status != MOI.CALLBACK_NODE_STATUS_INTEGER) && (return)
            eta_val = JuMP.callback_value(cb_data, eta)
            p_val = JuMP.callback_value(cb_data, p)
            y_val = JuMP.callback_value(cb_data, y)
            if (eta_val > 0.0) && (y_val > p_val + log(eta_val))
                inv_eta_val = 1/eta_val
                con = @build_constraint(y <= p + log(eta_val) + inv_eta_val*(eta - eta_val))
                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            end
        end
    
        get_ub(model) = JuMP.objective_value(model)
    
        MOI.set(model, MOI.LazyConstraintCallback(), outer_approximate)
    
        JuMP.optimize!(model)
        ub = get_ub(model)
        iteration = 0
    
        while (abs(ub - lb) > 0.01 || abs(ub - lb)/(1E-6 + abs(lb)) > 0.01)
            iteration += 1
            current_x_bus = Dict(i => JuMP.value(x_bus[i]) for i in interdictable_components.buses)
            current_x_line = Dict(i => JuMP.value(x_line[i]) for i in interdictable_components.lines)
            current_x_gen = Dict(i => JuMP.value(x_gen[i]) for i in interdictable_components.gens)
            current_buses = filter!(z -> last(z) > TOL, current_x_bus) |> keys |> collect
            current_lines = filter!(z -> last(z) > TOL, current_x_line) |> keys |> collect
            current_gens = filter!(z -> last(z) > TOL, current_x_gen) |> keys |> collect
            p_val = JuMP.value(p)
    
            cut_info = get_solution_values(input_data, current_gens, current_buses, current_lines)
            gp = map(y -> get(gen_prob, y, 1.0), current_gens)
            lp = map(y -> get(line_prob, y, 1.0), current_lines)
            bp = map(y -> get(bus_prob, y, 1.0), current_buses)
            probability = prod(gp) * prod(lp) * prod(bp)
            solution = Solution(current_lines, current_buses, current_gens, cut_info.load_shed * probability, cut_info.load_shed, probability, Dict())
            add_to_solutions(solutions, solution, max_num_solutions = solver_params["max_num_solutions"])
            
            # update lower bound
            sub_objective = cut_info.load_shed 
            sub_lb = p_val + log(sub_objective)
            (sub_lb > lb) && (lb = sub_lb) 
            
            # add solution removal cut 
            @constraint(model, sum([x_bus[i] for i in current_buses]; init = 0.0) + 
                sum([x_line[i] for i in current_lines]; init = 0.0) + 
                sum([x_gen[i] for i in current_gens]; init = 0.0) <= input_params["total_budget"] - 1)
            
            # add cutting plane 
            @constraint(model, eta <= sub_objective + 
                sum([cut_info.pg[i] * x_gen[i] for i in keys(cut_info.pg)]; init = 0.0) + 
                sum([cut_info.p[i] * x_line[i] for i in keys(cut_info.p)]; init = 0.0) +
                sum([cut_info.pb[i] * x_bus[i] for i in keys(cut_info.pb)]; init = 0.0)
            )
            
            MOI.set(model, MOI.LazyConstraintCallback(), outer_approximate)
            JuMP.optimize!(model)
            (get_ub(model) < ub) && (ub = get_ub(model))
            rel_gap = abs(ub - lb)/(1E-6 + abs(lb))
            print_table_cp(solutions |> length, iteration, 
                    get_best(solutions).cost, get_worst(solutions).cost, 
                    start_time, lb, ub, rel_gap, fields, field_chars)
            
            if time_limit_reached(start_time, solver_params["time_limit_in_seconds"]) || rel_gap <= 0.02
                break
            end 
    
        end
    
        return solutions
    end 