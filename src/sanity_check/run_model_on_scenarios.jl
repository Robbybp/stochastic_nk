using ArgParse
using HiGHS
using PowerModels
using JSON
using PrettyTables
using JuMP
using Gurobi
using Logging

# New ConsoleLogger that prints to stderr and accept messages with level >= Logging.Debug
debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger);

PMs = PowerModels
PMs.silence()

"""
julia --project=. src/sanity_check/run_model_on_scenarios.jl --case RTS_GMLC.m --scenario_file RTS_GMLC_1.json
"""

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--case", "-c"
        help = "case file"
        arg_type = String
        default = "pglib_opf_case240_pserc.m"

        "--scenario_file", "-s"
        help = "scenario file"
        arg_type = String
        default = "pglib_opf_case240_pserc_1.json"

        "--data_path", "-p"
        help = "data directory path"
        arg_type = String
        default = chop(Base.active_project(), tail = length("Project.toml")) * "data/"

        "--output_path"
        help = "output directory path"
        arg_type = String 
        default = chop(Base.active_project(), tail = length("Project.toml")) * "output/"
    end

    return parse_args(s)
end

function validate_parameters(params)
    mkpath(params["data_path"])
    mkpath(params["output_path"])
    case_file = params["data_path"] * "matpower/" * params["case"]
    if isfile(case_file) == false
        @error "$case_file does not exist, quitting."
        exit() 
    end  
    scenario_file = params["data_path"] * "scenario_data/" * params["scenario_file"]
    if isfile(scenario_file) == false
        @error "$scenario_file does not exist, quitting."
        exit() 
    end  
end 

function get_filenames_with_paths(params)
    matpower_file = params["data_path"] * "matpower/" * params["case"]
    scenario_file = params["data_path"] * "scenario_data/" * params["scenario_file"]
    return (mp_file = matpower_file, scenario_file = scenario_file)
end 

function main()
    cliargs = parse_commandline()
    # print the input parameters 
    pretty_table(cliargs, 
        title = "CLI parameters", 
        title_alignment = :c, 
        title_same_width_as_table = true, 
        show_header = false)

    validate_parameters(cliargs)
    files = get_filenames_with_paths(cliargs)
    run(cliargs, files)
    return 
end 

function run(cliargs, files)
    @info "starting run"
    mp_file = files.mp_file 
    scenario_file = files.scenario_file
    data = PowerModels.parse_file(mp_file; validate=false)
    PMs.make_per_unit!(data)
    data["total_load"] = [load["pd"] for (_, load) in data["load"]] |> sum
    for (_, gen) in get(data, "gen", [])
        (gen["pmax"] > 0.0) && (gen["pmin"] = 0.0)
        (gen["qmax"] > 0.0) && (gen["qmin"] = 0.0)
    end 
    ref = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    loads = ref[:load]
    scenario_data = JSON.parsefile(scenario_file)
    scenario_load_shed = Dict{String,Float64}(i => 0.0 for i in keys(scenario_data))

    for (id, scenario) in collect(scenario_data)
        @info "running scenario $id.."
        generators = scenario["gen"]
        lines = scenario["branch"]

        case = deepcopy(data)

        for i in generators 
            case["gen"][string(i)]["gen_status"] = 0
        end 
        for i in lines 
            case["branch"][string(i)]["br_status"] = 0
        end 
        PMs.propagate_topology_status!(case)
        load_shed, isolated_load_shed = run_dc_ls(case, loads)

        total_load_shed = isolated_load_shed + sum(values(load_shed))
        scenario_load_shed[id] = total_load_shed
        @debug "load shed = $total_load_shed"
    end 
    
    output_file = cliargs["output_path"] * "ls_" * cliargs["scenario_file"]
    
    open(output_file, "w") do f
        JSON.print(f, scenario_load_shed, 2)
    end
    return
end 

function run_dc_ls(case, loads)
    PMs.standardize_cost_terms!(case, order=2)
    PMs.calc_thermal_limits!(case)
    ref = PMs.build_ref(case)[:it][:pm][:nw][0]
    lp_optimizer = JuMP.optimizer_with_attributes(
            () -> HiGHS.Optimizer(), 
            "presolve" => "on",
            "output_flag" => false
        )
    lp_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer)
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

    @objective(model, Min, sum((1 - xd[i]) * load["pd"] for (i, load) in ref[:load]))
    
    @debug "number of ref bus constraints: $(length(ref[:ref_buses]))"
    for (i, _) in ref[:ref_buses]
        @constraint(model, va[i] == 0)
    end

    @debug "number of nodal balance constraints: $(length(ref[:bus]))"
    @debug ref[:bus_gens]
    for (i, _) in ref[:bus]
        # Build a list of the loads and shunt elements connected to the bus i
        bus_loads = [ref[:load][l] for l in ref[:bus_loads][i]]
        bus_shunts = [ref[:shunt][s] for s in ref[:bus_shunts][i]]

        # Active power balance at node i
        @constraint(model,
            sum(p_expr[a] for a in ref[:bus_arcs][i]) == # +                  
            sum(pg[g] for g in ref[:bus_gens][i]) -                 
            sum(xd[load["index"]] * load["pd"] for load in bus_loads) -                 
            sum(shunt["gs"] for shunt in bus_shunts)*1.0^2, 
            base_name="c_e_flow_equality_constr($i)"
        )
    end

    @debug "number of branch constraints: $(length(ref[:branch]))"
    for (i, branch) in ref[:branch]
        f_idx = (i, branch["f_bus"], branch["t_bus"])
        p_fr = p[f_idx]                     
        va_fr = va[branch["f_bus"]]         
        va_to = va[branch["t_bus"]]        
        _, b = PowerModels.calc_branch_y(branch)
        @constraint(model, branch["br_x"] * p_fr == (va_fr - va_to), base_name="c_e_flow_phase_constr($i)")
        @constraint(model, va_fr - va_to <= branch["angmax"], base_name="c_l_phase_diff_max_constr($i)")
        @constraint(model, va_fr - va_to >= branch["angmin"], base_name="c_l_phase_diff_min_constr($i)")
    end

    optimize!(model)
    
    # get the load shed for each individual load based on the solution
    existing_loads = ref[:load] |> keys 
    x_val = JuMP.value.(variables[:xd])
    all_loads = loads |> keys 
    load_shed = Dict(i => 0.0 for i in all_loads)
    isolated_load_shed = 0.0
    for i in all_loads 
        if !(i in existing_loads)
            isolated_load_shed +=  loads[i]["pd"] * 100.0
            continue 
        end 
        load_shed[i] = round((1-x_val[i]) * loads[i]["pd"] * 100.0; digits=4)
    end 
    @show load_shed, isolated_load_shed
    return load_shed, isolated_load_shed
end

main()