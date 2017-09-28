
function fetch_scenarios(config)
    scenariofile = config["scenariofile"]
    numscenarios = config["numscenarios"]
    numbatches = config["numbatches"]
    batchsize = config["batchsize"]
    batchid = config["batchid"]

    scenarios = readdlm(scenariofile)
    (size(scenarios)[1] < numscenarios) && (println(">> not enough scenarios available for the run, quitting program ...."); quit())
    (batchid > numbatches) && (println(">> number of batches and batch id inconsistent, quitting program ...."); quit())

    from = (batchid-1) * batchsize + 1
    to = batchid * batchsize
    
    return scenarios[from:to, :]
end
