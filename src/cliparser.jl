using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--case", "-c"
        help = "case file"
        arg_type = String
        default = "pglib_opf_case14_ieee.m"

        "--data_path", "-p"
        help = "data directory path"
        arg_type = String
        default = chop(Base.active_project(), tail = length("Project.toml")) * "data/"

        "--output_path"
        help = "output directory path"
        arg_type = String 
        default = chop(Base.active_project(), tail = length("Project.toml")) * "output/"

        "--problem", "-a"
        help = "problem selection - deterministic/stochastic"
        arg_type = String
        default = "deterministic"

        "--timeout", "-t"
        help = "time limit for the run in seconds"
        arg_type = Int
        default = 86400

        "--optimality_gap", "-o"
        help = "relative optimality gap in % (termination criteria)"
        arg_type = Float64
        default = 1e-2

        # the following options are valid only if the problem type is stochastic
        "--batch_id", "-b"
        help = "batch id for running the (id) batch"
        arg_type = Int
        default = 1

        "--num_batches", "-n"
        help = "total number of batches (num_batches should divide the total number of scenarios exactly)"
        arg_type = Int
        default = 10

        "--interdictable_components", "-i"
        help = "interdict lines or (lines and generators) (l/lg)"
        arg_type = String
        default = "l"

        "--budget", "-k"
        help = "budget for interdiction"
        arg_type = Int
        default = 2

        "--line_budget", "-l"
        help = "budget for lines"
        arg_type = Int 
        default = 2 

        "--generator_budget", "-g"
        help = "budget for generators"
        arg_type = Int 
        default = 0

        "--parallelize", "-x"
        help = "parallelize subproblem solves y/n"
        arg_type = String
        default = "n"

        "--workers", "-w"
        help = "number of workers"
        arg_type = Int
        default = 0
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
    budget_consistency = params["budget"] == params["generator_budget"] + params["line_budget"] 
    if budget_consistency == false 
        k = params["budget"] 
        g = params["generator_budget"] 
        l = params["line_budget"] 
        @error "line budget ($l) + generator budget ($g) does not equal the budget ($k), quitting."
        exit()
    end 
    # TODO: check for scenario and batch size consistency later for the stochastic interdiction
end 


function get_filenames_with_paths(params)
    case_name = chop(params["case"], tail = length(".m"))
    matpower_file = params["data_path"] * "matpower/" * params["case"]
    if params["problem"] == "deterministic"
        scenario_file = nothing
    else 
        scenario_file = params["data_path"] * "scenario_data/" * case_name * ".tar.gz"
        @info scenario_file
        if isfile(scenario_file) == false 
            @error "$scenario_file does not exist, quitting."
            exit()
        end 
    end 
    return (mp_file = matpower_file, scenario_file = scenario_file)
end 