using JSON 
using PowerModels 
using JuMP
using CPLEX

const PMs = PowerModels

PMs.silence()

function get_expected_load_sheds(file::String, 
    lines::Vector, 
    gens::Vector, 
    scenarios::Dict)::Dict 

    data = PowerModels.parse_file(file; validate=false)
    PowerModels.make_per_unit!(data)
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    bus_load = Dict{Int,Float64}()
    expected_bus_shed = Dict{Int,Float64}()
    bus_ratio = Dict{Int,Float64}()
    for (i, _) in ref[:bus]
        bus_loads = ref[:bus_loads][i]
        bus_shunts = ref[:bus_shunts][i]
        total_bus_load = sum(map(x -> ref[:load][x]["pd"], bus_loads); init=0.0)
        total_bus_shunt = sum(map(x -> ref[:shunt][x]["gs"], bus_shunts); init=0.0)
        bus_load[i] = total_bus_load + total_bus_shunt
        expected_bus_shed[i] = 0.0
    end 

    for (_, scenario) in scenarios 
        bus_shed = Dict{Int,Float64}()
        scenario_lines = scenario["branch"]
        scenario_gens = scenario["gen"]
        off_lines = unique!([scenario_lines..., lines...])
        off_gens = unique!([scenario_gens..., gens...])
        load_shed, shunt_shed = get_inner_solution(
            data, ref, off_gens, off_lines
        )
        for (i, _) in ref[:bus]
            bus_loads = ref[:bus_loads][i]
            bus_shunts = ref[:bus_shunts][i]
            total_bus_load_shed = sum(map(x -> load_shed[x], bus_loads); init=0.0)
            total_bus_shunt_shed = sum(map(x -> shunt_shed[x], bus_shunts); init=0.0)
            bus_shed[i] = total_bus_load_shed + total_bus_shunt_shed
            expected_bus_shed[i] += bus_shed[i]
        end 
    end 

    for (i, shed) in expected_bus_shed 
        expected_bus_shed[i] = shed/length(scenarios)
        bus_ratio[i] = (bus_load[i] == 0.0) ? 0.0 : expected_bus_shed[i]/bus_load[i]
    end 

    return bus_ratio
end 

""" Get load shed and power flow solution on interdictable components""" 
function get_inner_solution(data, ref, generators::Vector, lines::Vector; solver="cplex")::NamedTuple
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

    return run_dc_ls(case, ref)
end 

function run_dc_ls(case::Dict, original_ref::Dict; add_dc_lines_model::Bool=false)::NamedTuple 
    PowerModels.standardize_cost_terms!(case, order=2)
    PowerModels.calc_thermal_limits!(case)
    ref = PowerModels.build_ref(case)[:it][:pm][:nw][0]
    lp_optimizer = JuMP.optimizer_with_attributes(
        () -> CPLEX.Optimizer(), "CPX_PARAM_SCRIND" => 0
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
    @variable(model, 0 <= xs[i in keys(ref[:shunt])] <= 1)
    variables = Dict{Symbol,Any}(
        :va => va, 
        :pg => pg,
        :p => p,
        :xd => xd, 
        :xs => xs
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

    @objective(model, Min, 
        sum((1 - xd[i]) * load["pd"] for (i, load) in ref[:load]) +
        sum((1 - xs[i]) * shunt["gs"] for (i, shunt) in ref[:shunt]; init=0.0)
    )
    
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
                sum(xs[shunt["index"]] * shunt["gs"] for shunt in bus_shunts)*1.0^2
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
    existing_shunts = ref[:shunt] |> keys
    xd_val = JuMP.value.(variables[:xd])
    xs_val = JuMP.value.(variables[:xs])
    loads = original_ref[:load]
    shunts = original_ref[:shunt]
    all_loads = loads |> keys 
    all_shunts = shunts |> keys
    load_shed = Dict(i => 0.0 for i in all_loads)
    shunt_shed = Dict(i => 0.0 for i in all_shunts)
    for i in all_loads 
        if !(i in existing_loads)
            load_shed[i] = loads[i]["pd"] 
        else
            load_shed[i] = (1-xd_val[i]) * loads[i]["pd"]
        end 
    end 
    for i in all_shunts 
        if !(i in existing_shunts)
            shunt_shed[i] = shunts[i]["pd"]
        else
            shunt_shed[i] = (1-xs_val[i]) * shunts[i]["gs"]
        end 
    end 

    return (load_shed = load_shed, shunt_shed = shunt_shed)
end 

function main()
    root_dir = chop(Base.active_project(), tail = length("Project.toml"))
    scenario_file = root_dir * "data/scenario_data/RTS_GMLC_1.json"
    case_file = root_dir * "data/matpower/RTS_GMLC.m"
    scenarios = JSON.parsefile(scenario_file)
    result_dir = root_dir * "output/stochastic/"

    to_write = Dict()
    for k in 1:10 
        result_file = result_dir * "RTS-GMLC--sto--RTS-GMLC-1--k" * string(k) * "--m200.json"
        results = JSON.parsefile(result_file)["results"]
        lines = results["lines"]
        generators = results["generators"]
        to_write[k] = get_expected_load_sheds(case_file, lines, generators, scenarios)
    end 
         
    file = root_dir * "/visualization/data/expected_ls_factor.json" 
    open(file, "w") do f 
        JSON.print(f, to_write, 2) 
    end
end 

main() 