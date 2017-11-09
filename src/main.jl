using JuMP
using PowerModels
using CPLEX
using Gurobi
using MathProgBase
 
PMs = PowerModels

include("parse.jl")
include("utils.jl")
include("full_model.jl")
include("bound_tightening.jl")
include("master.jl")
include("subproblem.jl")
include("pdkkt.jl")
include("benders.jl")

# setting up configuration for the run
config = parse_commandline()
config["casefile"] = string(config["path"], config["file"])
config["numscenarios"] = config["batchsize"] * config["numbatches"]
solver_to_use = CplexSolver(CPX_PARAM_THREADS=1)
(config["solver"] == "gurobi") && (solver_to_use = GurobiSolver(OutputFlag=0,Threads=1))

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

if config["parallel"] == "y"
    addprocs(config["workers"])
    @everywhere using MathProgBase
    @everywhere using CPLEX
end

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
    @time solve(m)
    println(">> obj full model: $(getobjectivevalue(m))")
    println(">> $(getvalue(getindex(m, :x)))")
    println(">> $(getvalue(getindex(m, :y)))")    
    
end

if config["algo"] != "full"
    
    if config["bt"] == "y"
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
    end
    
    # compute master problem auxilliary variable upper bounds
    println(">> computing tight bounds for master problem auxiliary variables")
    config["theta_ub"] = Dict{Int,Float64}()
    for s in 1:config["batchsize"]
        m = post_dc_primal(data, scenarios[s,:], Model(solver=CplexSolver(CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0)))
        solve(m)
        config["theta_ub"][s] = getobjectivevalue(m)
    end

    # create the relaxed master problem for the first iteration 
    println(">> creating master problem")
    master = create_master_model(scenarios, ref, config, Model(solver=solver_to_use))
    
    # create subproblem nonchanging matrices
    A, sense, l, u = create_matrices(scenarios, ref, config, Model())

    # create rhs variable expression vector for master cuts
    master.ext = create_expression_vectors(scenarios, ref, config, master_model=master)
    
    println(">> Lshaped")
    tic()
    @time status, obj, sol = Lshaped(scenarios, ref, config, A, sense, l, u, master, CplexSolver(CPX_PARAM_THREADS=1,CPX_PARAM_SCRIND=0))
    config["time"] = toq()
    println(">> algorithm ended")

    @printf ">> obj = %.4f MW \n" obj*ref[:nw][0][:baseMVA]
    println(">> branch_indexes = $(sol[:x])")
    println(">> gen_indexes = $(sol[:y])")
    
    config["obj"] = obj*ref[:nw][0][:baseMVA]
    config["sol"] = sol
    total_load = 0
    for (i, bus) in ref[:nw][0][:bus]
        total_load += bus["pd"]
    end
    total_load *= ref[:nw][0][:baseMVA]
    
    config["percent_load_shed"] = config["obj"]/total_load*100

    write_solution(config, ref)

end
