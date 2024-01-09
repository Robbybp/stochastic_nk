using PowerModels 
using JuMP 
using Gurobi 
using PrettyTables
using Dates
using JSON
using Logging
using UUIDs

const GRB_ENV = Gurobi.Env()

PowerModels.silence()

include("cliparser.jl")
include("data_helper.jl")
include("types.jl")
include("common/dc-ls.jl")
include("deterministic/run.jl")
include("stochastic/run.jl")

# New ConsoleLogger that prints to stderr and accept messages with level >= Logging.Debug
debug_logger = ConsoleLogger(stderr, Logging.Debug)
global_logger(debug_logger)

# timing functions 
time_limit_reached(start_time::DateTime, limit_in_seconds::Float64) = 
    ((now() - start_time).value/1000.0) |> round > limit_in_seconds

get_time(start_time) = ((now() - start_time).value/1000.0) |> round

function main()
    config = parse_commandline()
    
    # print the input parameters 
    pretty_table(config, 
        title = "CLI parameters", 
        title_alignment = :c, 
        title_same_width_as_table = true, 
        show_header = false)

    validate_parameters(config)
    files = get_filenames_with_paths(config)
    run(config, files)
    return 
end 

function run(config, files)
    mp_file = files.mp_file 
    scenario_file = files.scenario_file 
    
    if config["problem"] == "deterministic"
        results = run_deterministic(config, mp_file) 
        write_results(config, results)
    end
    
    if config["problem"] == "stochastic"
        results = run_stochastic(config, mp_file, scenario_file) 
        write_results(config, results)
    end
end 

function get_config_data_dict(config::Dict)
    config_data = Dict(
        "case" => config["case"],
        "problem" => config["problem"],
        "budget" => config["budget"],
        "separate_budgets" => config["use_separate_budgets"]
    )
    if config["use_separate_budgets"]
        config_data["line_budget"] = config["line_budget"]
        config_data["generator_budget"] = config["generator_budget"]
    else 
        config_data["line_budget"] = NaN
        config_data["generator_budget"] = NaN
    end 
    if config["problem"] == "stochastic"
        config_data["scenario_data"] = config["scenario_file"]
        config_data["num_scenarios"] = config["maximum_scenarios"]
    else 
        config_data["scenario_data"] = ""
        config_data["num_scenarios"] = NaN
    end 
    return config_data
end   

function get_run_data_dict(results::Results)
    return Dict(
        "time_ended" => string(now()), 
        "objective" => round(results.objective_value; digits=4), 
        "bound" => round(results.bound; digits=4), 
        "run_time" => round(results.run_time_in_seconds; digits=2), 
        "relative_gap" => round(results.optimality_gap; digits=2), 
        "lines" => results.solution.lines, 
        "generators" => results.solution.generators
    )
end 

function write_results(config::Dict, results::Results)
    config_data = get_config_data_dict(config)
    run_data = get_run_data_dict(results)

    to_write = Dict("instance_data" => config_data, "results" => run_data)
    file = config["output_path"] * config_data["problem"] * "/" * string(uuid1()) * ".json"
    open(file, "w") do f 
        JSON.print(f, to_write, 2)
    end 
end 

main()