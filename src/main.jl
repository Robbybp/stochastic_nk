using JuMP
using PowerModels
using CPLEX
 
PMs = PowerModels

include("parse.jl")
include("utils.jl")
include("full_model.jl")
include("subproblem.jl")

PMs = PowerModels
 
include("parse.jl")
include("utils.jl")
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

if (config["algo"] == "full")
    #mp = post_dc_primal(data, scenarios, Model(solver=CplexSolver())) 
    #solve(mp)
    #println(">> obj primal first scenario: $(getobjectivevalue(mp))")
    #md = post_dc_dual(data, scenarios, Model(solver=CplexSolver())) 
    #solve(md)
    #println(">> obj dual first scenario: $(getobjectivevalue(md))")
    kkt = post_dc_kkt(data, scenarios, Model(solver=CplexSolver())) 
    solve(kkt)
    println(">> obj kkt first scenario: $(getobjectivevalue(kkt))")
    m = create_full_model(scenarios, ref, config, model=Model(solver=CplexSolver()))
    solve(m)
    println(">> obj: $(getobjectivevalue(m))")
    println(">> $(getvalue(getindex(m, :x)))")
    println(">> $(getvalue(getindex(m, :y)))")
    
    println(kkt)
    println(scenarios)
    println(m)
end
