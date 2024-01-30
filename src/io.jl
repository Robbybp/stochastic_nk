""" generate config_data dictionary to write to file from cli-args """
function get_config_data(cliargs::Dict)
    config_data = Dict(
        "case" => cliargs["case"],
        "problem" => cliargs["problem"],
        "budget" => cliargs["budget"],
        "separate_budgets" => cliargs["use_separate_budgets"]
    )
    if cliargs["use_separate_budgets"]
        config_data["line_budget"] = cliargs["line_budget"]
        config_data["generator_budget"] = cliargs["generator_budget"]
    else 
        config_data["line_budget"] = NaN
        config_data["generator_budget"] = NaN
    end 
    if cliargs["problem"] in ["stochastic", "evp"]
        config_data["scenario_data"] = cliargs["scenario_file"]
        config_data["num_scenarios"] = cliargs["maximum_scenarios"]
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
    lgk = if cd["separate_budgets"]
        "lk" * string(cd["line_budget"]) * "gk" * string(cd["generator_budget"])
        else ""
        end 
    s = isempty(cd["scenario_data"]) ? "NaN" : replace(first(split(cd["scenario_data"], ".")), "_" => "-")
    m = "m" * string(cd["num_scenarios"])
    (p == "det") && (return join(filter(!=(""), [c, p, k, lgk]), "--") * ".json")
    return join(filter(!=(""), [c, p, s, k, lgk, m]), "--") * ".json"
end 

""" generate run info dictionary to write to file """
function get_run_data(results::Results)
    return Dict(
        "time_ended" => string(now()), 
        "objective" => round(results.objective_value * 100.0; digits=4), 
        "bound" => round(results.bound * 100.0; digits=4), 
        "run_time" => round(results.run_time_in_seconds; digits=2), 
        "relative_gap" => round(results.optimality_gap; digits=2), 
        "lines" => results.solution.lines, 
        "generators" => results.solution.generators
    )
end 

""" write results to file """
function write_results(cliargs::Dict, results::Results)
    config_data = get_config_data(cliargs)
    run_data = get_run_data(results)

    to_write = Dict("instance_data" => config_data, "results" => run_data)
    file = cliargs["output_path"] * config_data["problem"] * "/" * get_outfile_name(config_data)
    open(file, "w") do f 
        JSON.print(f, to_write, 2)
    end 
end 
