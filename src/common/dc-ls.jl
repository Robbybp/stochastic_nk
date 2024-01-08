""" Get load shed and power flow solution on interdictable components""" 
function get_inner_solution(data, ref, generators::Vector, lines::Vector; use_pm::Bool=false)::NamedTuple
    case_data = data
    # deepcopy and turn-off interdicted components 
    case = deepcopy(case_data)
    for i in generators 
        case["gen"][string(i)]["gen_status"] = 0
    end 
    for i in lines 
        case["branch"][string(i)]["br_status"] = 0
    end 
    PowerModels.propagate_topology_status!(case)

    if use_pm
        lp_optimizer = JuMP.optimizer_with_attributes(
            () -> Gurobi.Optimizer(GRB_ENV), "LogToConsole" => 0
        )
        pm = instantiate_model(case, DCPPowerModel, PowerModels._build_mld)
        result = optimize_model!(pm, optimizer = lp_optimizer)

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
    return run_dc_ls(case, ref)
end 

function run_dc_ls(case::Dict, original_ref::Dict; add_dc_lines_model::Bool=false)::NamedTuple 
    PowerModels.standardize_cost_terms!(case, order=2)
    PowerModels.calc_thermal_limits!(case)
    ref = PowerModels.build_ref(case)[:it][:pm][:nw][0]
    lp_optimizer = JuMP.optimizer_with_attributes(
        () -> Gurobi.Optimizer(GRB_ENV), "LogToConsole" => 0
    )
    model = Model(lp_optimizer)
    
    @variable(model, va[i in keys(ref[:bus])])
    @variable(model, 
        ref[:gen][i]["pmin"] <= 
        pg[i in keys(ref[:gen])] <= 
        ref[:gen][i]["pmax"]
    )
    @variable(model, 
        -ref[:branch][l]["rate_a"] <= 
        p[(l,i,j) in ref[:arcs_from]] <= 
        ref[:branch][l]["rate_a"]
    )
    p_expr = Dict([((l,i,j), 1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]])
    p_expr = merge(p_expr, Dict([((l,j,i), -1.0*p[(l,i,j)]) for (l,i,j) in ref[:arcs_from]]))
    @variable(model, 0 <= xd[i in keys(ref[:load])] <= 1)
    variables = Dict{Symbol,Any}(
        :va => va, 
        :pg => pg,
        :p => p,
        :xd => xd
    ) 
    if (add_dc_lines_model)
        @variable(model, p_dc[a in ref[:arcs_dc]])
        variables[:p_dc] = p_dc
        for (l,dcline) in ref[:dcline]
            f_idx = (l, dcline["f_bus"], dcline["t_bus"])
            t_idx = (l, dcline["t_bus"], dcline["f_bus"])
    
            JuMP.set_lower_bound(p_dc[f_idx], dcline["pminf"])
            JuMP.set_upper_bound(p_dc[f_idx], dcline["pmaxf"])
    
            JuMP.set_lower_bound(p_dc[t_idx], dcline["pmint"])
            JuMP.set_upper_bound(p_dc[t_idx], dcline["pmaxt"])
        end
        for (i, dcline) in ref[:dcline]
            f_idx = (i, dcline["f_bus"], dcline["t_bus"])
            t_idx = (i, dcline["t_bus"], dcline["f_bus"])  
            @constraint(model, 
                (1-dcline["loss1"])*p_dc[f_idx] + (p_dc[t_idx] - dcline["loss0"]) == 0,
                base_name = "c_dc_line($i)"
                )
        end
    end 

    @objective(model, Min, sum((1 - xd[i]) * load["pd"] for (i, load) in ref[:load]))
    
    for (i, _) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    for (i, _) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        if (add_dc_lines_model)
            p_dc = get(variables, :p_dc, nothing)
            @constraint(model,
                sum(p_expr[a] for a in ref[:bus_arcs][i]) +                  
                sum(p_dc[a_dc] for a_dc in ref[:bus_arcs_dc][i]) ==     
                sum(pg[g] for g in ref[:bus_gens][i]) -                 
                sum(xd[load["index"]] * load["pd"] for load in bus_loads) -                 
                sum(shunt["gs"] for shunt in bus_shunts)*1.0^2
            )
            continue 
        end 

        # Active power balance at node i
        @constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) ==                   
            sum(pg[g] for g in ref[:bus_gens][i]) -                 
            sum(xd[load["index"]] * load["pd"] for load in bus_loads) -                 
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2, 
            base_name="c_e_flow_equality_constr($i)"
        )
    end

    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        p_fr = p[f_idx]                     
        va_fr = va[branch["f_bus"]]         
        va_to = va[branch["t_bus"]]        
        # _, _ = PowerModels.calc_branch_y(branch)
        @constraint(model, branch["br_x"] * p_fr == (va_fr - va_to), base_name="c_e_flow_phase_constr($i)")
        @constraint(model, va_fr - va_to <= branch["angmax"], base_name="c_l_phase_diff_max_constr($i)")
        @constraint(model, va_fr - va_to >= branch["angmin"], base_name="c_l_phase_diff_min_constr($i)")
    end
    optimize!(model)

    # get the load shed for each individual load based on the solution
    existing_loads = ref[:load] |> keys 
    x_val = JuMP.value.(variables[:xd])
    loads = original_ref[:load]
    all_loads = loads |> keys 
    load_shed = Dict(i => 0.0 for i in all_loads)
    isolated_load_shed = 0.0
    for i in all_loads 
        if !(i in existing_loads)
            isolated_load_shed +=  loads[i]["pd"]
            continue 
        end 
        load_shed[i] = (1-x_val[i]) * loads[i]["pd"]
    end 
    total_load_shed = isolated_load_shed + sum(values(load_shed))
    pg_values = Dict(i => JuMP.value(pg[i]) for i in keys(ref[:gen]))
    p_values = Dict(l => abs(JuMP.value(p[(l, i, j)])) for (l, i, j) in ref[:arcs_from])

    return (load_shed = total_load_shed, pg = pg_values, p = p_values)
end 


""" get load shed and power flow solution on fractional interdictable components""" 
# function get_inner_solution(data, ref, generators::Dict{Int,Float64}, lines::Dict{Int,Float64})::NamedTuple
#     case_data = data
#     case = deepcopy(case_data)
#     for i in generators 
#         case["gen"][string(i)]["gen_status"] = 0
#     end 

#     for i in lines 
#         case["branch"][string(i)]["br_status"] = 0
#     end 

#     PowerModels.propagate_topology_status!(case)
#     lp_optimizer = JuMP.optimizer_with_attributes(
#         () -> Gurobi.Optimizer(GRB_ENV), "LogToConsole" => 0
#     )
    
#     pm = instantiate_model(case, DCPPowerModel, PowerModels._build_mld)
#     result = optimize_model!(pm, optimizer = lp_optimizer)

#     # result = solve_model(case, DCPPowerModel, lp_optimizer, PowerModels._build_mld)
#     load_served = [load["pd"] for (_, load) in result["solution"]["load"]] |> sum
#     load_shed = case_data["total_load"] - load_served

#     pg = Dict(i => result["solution"]["gen"][string(i)]["pg"]
#         for i in keys(ref[:gen]) if haskey(result["solution"]["gen"], string(i)))
    
#     p = Dict(i => max(
#         abs(result["solution"]["branch"][string(i)]["pf"]), 
#         abs(result["solution"]["branch"][string(i)]["pt"])
#         ) for i in keys(ref[:branch]) if haskey(result["solution"]["branch"], string(i)) 
#     )

#     return (load_shed = load_shed, pg = pg, p = p)
# end 
