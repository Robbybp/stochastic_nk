using Distributions
using PowerModels
using ArgParse

PMs = PowerModels


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--ifile", "-i"
        help = "pglib case file name with additional fields for n-k"
        arg_type = String
        default = "pglib_opf_case14_ieee_nk.m"
        
        "--ipath", "-p"
        help = "data directory path"
        arg_type = String
        default = "../data/"

        "--numscenarios", "-n"
        help = "number of scenarios to be generated"
        arg_type = Int
        default = 50000
        
        "--ofile", "-o"
        help = "output file name"
        arg_type = String
        default = "case14_scenarios.txt"

        "--opath", "-a"
        help = "output file path"
        arg_type = String
        default = "../data/scenario_data/"

    end

    return parse_args(s)
end


config = parse_commandline()
config["inputfile"] = string(config["ipath"], config["ifile"])
config["outputfile"] = string(config["opath"], config["ofile"])

(isfile(config["outputfile"])) && (println(">> scenario file exists ... exiting program"); quit())

data = PMs.parse_file(config["inputfile"])
zipfilename = string(config["opath"], "case", length(data["bus"]), ".tar.gz")
(isfile(zipfilename)) && (println(">> scenario zipped file exists ... exiting program"); quit())

ref = PMs.build_ref(data)[:nw][0]
numcolumns = length(keys(ref[:gen])) + length(keys(ref[:branch]))
scenarios = zeros(config["numscenarios"], numcolumns)

dist = Dict( [ (:gen, Dict{Int,Distributions.Bernoulli{Float64}}()), (:branch, Dict{Int,Distributions.Bernoulli{Float64}}()) ] )
dist[:gen] = Dict( [ (i, Bernoulli(gen["prob"])) for (i, gen) in ref[:gen] ] )
dist[:branch] = Dict( [ (i, Bernoulli(branch["prob"])) for (i, branch) in ref[:branch] ] )

for i in 1:config["numscenarios"]
    scenario = Int[]
    for g in keys(ref[:gen])
        push!(scenario, rand(dist[:gen][g]))
    end
    for l in keys(ref[:branch])
        push!(scenario, rand(dist[:branch][l]))
    end
    @assert length(scenario) == numcolumns
    scenarios[i,:] = scenario
end

scenarios = round.(Int, scenarios)

writedlm(config["outputfile"], scenarios, " ")
run(`tar -zcvf $zipfilename $(config["outputfile"])`)
run(`rm -f $(config["outputfile"])`)
