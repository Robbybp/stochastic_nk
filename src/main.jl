using PowerModels 
using JuMP 
using Gurobi 
using PrettyTables
using Dates
using JSON
using Logging

const GRB_ENV = Gurobi.Env()

PowerModels.silence()

include("cliparser.jl")
include("data_helper.jl")
include("common/dc-ls.jl")
include("deterministic/types.jl")
include("deterministic/run.jl")
include("stochastic/types.jl")
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
        @show results
        # TODO: implement result writing to file 
        # write_results(config, results)
    end
    # TODO: need to implement 
    
    # if config["problem"] == "stochastic"
    #     results = run_stochastic(config, mp_file, scenario_file) 
    #     # write_results(config, results)
    # end
end 

main()