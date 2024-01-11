""" generate config dictionary to write to file from cli-args """
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

""" generate output file name from the config_data dictionary """
function get_outfile_name(cd::Dict)
    c = replace(first(split(cd["case"], ".")), "_" => "-")
    p = cd["problem"][1:3]
    k = "k" * string(cd["budget"])
    sk = "sk" * string(cd["separate_budgets"])
    lk = "lk" * string(cd["line_budget"]) 
    gk = "gk" * string(cd["generator_budget"])
    s = isempty(cd["scenario_data"]) ? "NaN" : replace(first(split(cd["scenario_data"], ".")), "_" => "-")
    m = "m" * string(cd["num_scenarios"])
    return join([c, p, s, k, sk, lk, gk, m, ".json"], "--")
end 

""" generate run info dictionary to write to file """
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

""" write results to file """
function write_results(config::Dict, results::Results)
    config_data = get_config_data_dict(config)
    run_data = get_run_data_dict(results)

    to_write = Dict("instance_data" => config_data, "results" => run_data)
    file = config["output_path"] * config_data["problem"] * "/" * get_outfile_name(config_data)
    open(file, "w") do f 
        JSON.print(f, to_write, 2)
    end 
end 
