using JuMP
using PowerModels
using CPLEX
using Gurobi

 
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
solver_to_use = CplexSolver()
(config["solver"] == "gurobi") && (solver_to_use = GurobiSolver(OutputFlag=0))

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
        mp = post_dc_primal(data, scenarios, Model(solver=solver_to_use)) 
        solve(mp)
        println(">> obj primal one scenario: $(getobjectivevalue(mp))")

        data = PMs.parse_file(config["casefile"])
        md = post_dc_dual(data, scenarios, Model(solver=solver_to_use)) 
        solve(md)
        println(">> obj dual one scenario: $(getobjectivevalue(mp))")
    
        data = PMs.parse_file(config["casefile"])
        mkkt = post_dc_kkt(data, scenarios, Model(solver=solver_to_use)) 
        solve(mkkt)
        println(">> obj kkt one scenario: $(getvalue(mkkt.ext[:primalobj_expr]))")
        
    end

    data = PMs.parse_file(config["casefile"])
    m = create_full_model(scenarios, ref, config, Model(solver=solver_to_use))
    solve(m)
    println(">> obj full model: $(getobjectivevalue(m))")
    println(">> $(getvalue(getindex(m, :x)))")
    println(">> $(getvalue(getindex(m, :y)))")    
    
end

if config["algo"] == "Lshaped" || config["algo"] == "Lshapedreg"
    println(">> creating and solving full LP relaxation")
    # solve relaxation of the full model (to be used for bound-tightening)
    full_model = create_full_model(scenarios, ref, config, Model(solver=solver_to_use))
    solve(full_model, relaxation=true)
    relaxation_obj = getobjectivevalue(full_model)
    println(">> relaxation solved, objective value: $relaxation_obj")
    
    # perform bound tightening 
    println(">> tightening bounds")
    bt_model = create_bound_tightening_model(scenarios, ref, config, Model(solver=solver_to_use), relaxation_obj=relaxation_obj)
    tighten_bounds(scenarios, ref, config, bt_model)
    println(">> bounds tightened")
    
    # resolve relaxation with tightened bounds as a check
    println(">> resolving relaxation with tightened bounds")
    full_model = create_full_model(scenarios, ref, config, Model(solver=solver_to_use))
    solve(full_model, relaxation=true)
    println(">> resolved objective value: $(getobjectivevalue(full_model))")
    @assert isapprox(relaxation_obj, getobjectivevalue(full_model), atol=1e-6)
    println(">> check for correctness of bound-tightening cleared")


    # create the relaxed master problem for the first iteration 
    println(">> creating master problem")
    master = create_master_model(scenarios, ref, config, Model(solver=solver_to_use))   
    

end
