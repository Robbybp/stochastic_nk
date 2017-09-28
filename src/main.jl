 using JuMP
 using PowerModels
 using CPLEX
 
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
 
 data = PMs.parse_file("../data/pglib_opf_case24_ieee_rts.m")
 m_primal = post_dc_primal(data, Model(solver=CplexSolver()))
 
 status = solve(m_primal)
 
 m_dual = post_dc_dual(data, Model(solver=CplexSolver()))
 
 status = solve(m_dual)
 
 @assert getobjectivevalue(m_primal) == getobjectivevalue(m_dual)
