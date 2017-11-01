using JuMP
using PowerModels
using CPLEX
 
PMs = PowerModels

include("parse.jl")
include("utils.jl")
include("full_model.jl")
include("bound_tightening.jl")
include("master.jl")
include("subproblem.jl")

# setting up configuration for the run
config = parse_commandline()
config["casefile"] = string(config["path"], config["file"])
config["numscenarios"] = config["batchsize"] * config["numbatches"]

for (arg, val) in config
    println(">> $arg => $val") 
end

data = PMs.parse_file(config["casefile"])
ref = PMs.build_ref(data)

config["scenariofile"] = string(config["path"], "scenario_data/case", length(keys(ref[:nw][0][:bus])), "_scenarios_unique.txt")

# fetch the scenarios for the run
(!isfile(config["scenariofile"])) && (println(">> scenario file does not exist, quitting program ..."); quit())
scenarios = fetch_scenarios(config)
@assert size(scenarios)[1] == config["batchsize"]

if config["algo"] == "full"
    if config["batchsize"] == 1
        mp = post_dc_primal(data, scenarios, Model(solver=CplexSolver())) 
        solve(mp)
        println(">> obj primal one scenario: $(getobjectivevalue(mp))")
    end

    m = create_full_model(scenarios, ref, config, Model(solver=CplexSolver()))
    solve(m)
    println(">> obj full model: $(getobjectivevalue(m))")
    println(">> $(getvalue(getindex(m, :x)))")
    println(">> $(getvalue(getindex(m, :y)))")    
end

if config["algo"] == "Lshaped" || config["algo"] == "Lshapedreg"
    
    # perform bound tightening 
    bt_model = create_bound_tightening_model(scenarios, ref, config, Model(solver=CplexSolver()))
    tighten_bounds(scenarios, ref, config, bt_model)
    println(scenarios)
    println(config["bounds"][:primal_dclb][1])
    println(config["bounds"][:primal_dcub][1])

    # create the relaxed master problem for the first iteration 
    master = create_master_model(scenarios, ref, config, Model(solver=CplexSolver()))   
    

end
