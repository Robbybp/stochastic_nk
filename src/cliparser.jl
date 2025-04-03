using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--case", "-c"
        help = "case file"
        arg_type = String
        default = "pglib_opf_case240_pserc.m"

        "--data_path", "-p"
        help = "data directory path"
        arg_type = String
        default = chop(Base.active_project(), tail = length("Project.toml")) * "data/"

        "--output_path"
        help = "output directory path"
        arg_type = String 
        default = chop(Base.active_project(), tail = length("Project.toml")) * "output/"

        "--problem", "-a"
        help = "problem selection - deterministic/stochastic/evp"
        arg_type = String
        default = "stochastic"

        "--timeout", "-t"
        help = "time limit for the run in seconds"
        arg_type = Int
        default = 86400

        "--optimality_gap", "-o"
        help = "relative optimality gap in % (termination criteria)"
        arg_type = Float64
        default = 1e-2

        "--budget", "-k"
        help = "budget for interdiction"
        arg_type = Int
        default = 2

        "--use_separate_budgets" 
        help = "use separate line and generator budgets" 
        action = :store_true 

        "--bus_budget"
        help = "budget for buses. CANNOT BE USED WITH --line_budget OR --generator_budget"
        arg_type = Int 
        default = 0

        "--line_budget", "-l"
        help = "budget for lines"
        arg_type = Int 
        default = 2 

        "--generator_budget", "-g"
        help = "budget for generators"
        arg_type = Int 
        default = 0

        "--inner_solver"
        help = "cplex/gurobi"
        arg_type = String 
        default = "cplex"

        "--rerun"
        help = "re-run even if result file already exists"
        action = :store_true

        # the following options are valid only if the problem type is stochastic
        "--scenario_file", "-s"
        help = "scenario file"
        arg_type = String
        default = "pglib_opf_case240_pserc_1.json"

        "--maximum_scenarios", "-m"
        help = "limits the total number of scenarios used from scenario file"
        arg_type = Int
        default = 50
    end

    return parse_args(s)
end

function validate_parameters(params; skip_path_validation::Bool = false)
    # I want a version of this code that does not interact so much with the
    # filesystem. I will handle IO myself at a higher level.
    # TODO: How exactly should I handle this slightly different interface?
    if skip_path_validation
        mkpath(params["data_path"])
        mkpath(params["output_path"])
        case_file = params["data_path"] * "matpower/" * params["case"]
        if isfile(case_file) == false
            @error "$case_file does not exist, quitting."
            exit() 
        end  
    end
    if params["bus_budget"] > 0 && params["generator_budget"] + params["line_budget"] > 1
        error("Cannot specify a bus budget and a generator or line budget")
    end
    if (params["use_separate_budgets"])
        if params["bus_budget"] > 0
            error("Cannot specify a bus budget while using separate budgets")
        end
        budget_consistency = params["total_budget"] == params["generator_budget"] + params["line_budget"] + params["bus_budget"]
        if budget_consistency == false 
            k = params["budget"] 
            b = params["bus_budget"] 
            g = params["generator_budget"] 
            l = params["line_budget"] 
            @error "line budget ($l) + generator budget ($g) + bus budget ($b) does not equal the budget ($k), quitting."
            exit()
        end 
    end 
    # TODO: check for scenario and batch size consistency later for the stochastic interdiction
end 


function get_filenames_with_paths(params)
    case_name = chop(params["case"], tail = length(".m"))
    matpower_file = params["data_path"] * "matpower/" * params["case"]
    if params["problem"] == "deterministic"
        scenario_file = nothing
    else 
        scenario_file = params["data_path"] * "scenario_data/" * params["scenario_file"]
    end 
    return (mp_file = matpower_file, scenario_file = scenario_file)
end 
