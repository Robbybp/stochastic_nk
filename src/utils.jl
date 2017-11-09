
function fetch_scenarios(config, numbuses)
    scenariofile = config["scenariofile"]
    zipfile = config["zipfile"]
    numscenarios = config["numscenarios"]
    numbatches = config["numbatches"]
    batchsize = config["batchsize"]
    batchid = config["batchid"]
    scenario_dir = string(config["path"], "scenario_data/")
    
    (!isfile(scenariofile)) && (run(`tar -zxvf $zipfile -C $scenario_dir`))
    scenarios = readdlm(scenariofile)
    (isfile(zipfile)) && (run(`rm -f $scenariofile`))

    (size(scenarios)[1] < numscenarios) && (println(">> not enough scenarios available for the run, quitting program ...."); quit())
    (batchid > numbatches) && (println(">> number of batches and batch id inconsistent, quitting program ...."); quit())

    from = (batchid-1) * batchsize + 1
    to = batchid * batchsize
    
    return scenarios[from:to, :]
end

function write_solution(config, ref)
    println(">> writing solution to file")
    delete!(config, "bounds")
    delete!(config, "theta_ub")
    filename = string("../output_from_runs/case", length(ref[:nw][0][:bus]), "_s", config["batchsize"], "_id", config["batchid"], "_k", config["budget"], ".txt")
    println(">> output filename: $filename")

    open(filename, "w") do f
        for (i, val) in config
            if i == "sol"
                write(f, "branch = $(sol[:x])\n")
                write(f, "gen = $(sol[:y])\n")
            else
                write(f, "$i = $val\n")
            end
        end
    end

    return
end

