using Distributions
using PowerModels
using ArgParse
using PrettyTables
using DataFrames
using CSV
using Random
using JSON

PMs = PowerModels
PMs.silence()


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--case", "-c"
        help = "case file name"
        arg_type = String
        default = "pglib_opf_case14_ieee.m"
        
        "--data_path", "-p"
        help = "data directory path"
        arg_type = String
        default = chop(Base.active_project(), tail = length("Project.toml")) * "data/"

        "--output_path"
        help = "output directory path"
        arg_type = String 
        default = chop(Base.active_project(), tail = length("Project.toml")) * "data/scenario_data/"

        "--num_scenarios", "-n"
        help = "number of scenarios to be generated"
        arg_type = Int
        default = 200

        "--num_outages", "-l"
        help = "maximum number of line/gen outages in each scenario"
        arg_type = Int 
        default = 4 

        "--use_clusters"
        help = "flag to use the clustering data for the case file"
        action = :store_true
    end

    return parse_args(s)
end

function main()
    config = parse_commandline()
    
    # print the input parameters 
    pretty_table(config, 
        title = "CLI parameters", 
        title_alignment = :c, 
        title_same_width_as_table = true, 
        show_header = false)

    validate_parameters(config)
    files = get_filenames_with_paths(config)
    run(config, files)
    return 
end 

function validate_parameters(params)
    mkpath(params["data_path"])
    mkpath(params["output_path"])
    case_file = params["data_path"] * "matpower/" * params["case"]
    if isfile(case_file) == false
        @error "$case_file does not exist, quitting."
        exit() 
    end  
end 


function get_filenames_with_paths(params)
    case_name = chop(params["case"], tail = length(".m"))
    matpower_file = params["data_path"] * "matpower/" * params["case"]
    scenario_file = params["data_path"] * "scenario_data/" * case_name * ".json"
    tar_file = params["data_path"] * "scenario_data/" * case_name * ".tar.gz"
    cluster_file = params["data_path"] * "gis/" * case_name * "_cluster.csv"
    return (mp_file = matpower_file, 
        scenario_file = scenario_file, 
        cluster_file = cluster_file,
        tar_file = tar_file)
end 

function run(config, files)
    mp_file = files.mp_file 
    scenario_file = files.scenario_file 
    cluster_file = files.cluster_file
    tar_file = files.tar_file
    (isfile(scenario_file)) && (@error "$scenario_file exists, quitting"; exit())
    (config["use_clusters"] && !isfile(cluster_file)) && (@error "$cluster_file does not exists, quitting"; exit())

    data = PMs.parse_file(mp_file)
    ref = PMs.build_ref(data)[:it][:pm][:nw][0]
    @info "number of buses: $(length(ref[:bus]))"
    Random.seed!(0)

    if config["use_clusters"] == false
        line_ids = ref[:branch] |> keys 
        gen_ids = ref[:gen] |> keys
        components = ["line", "gen"]

        num_outages = config["num_outages"]
        scenarios = Dict()
        for i in 1:config["num_scenarios"]
            scenarios[i] = Dict("branch" => [] , "gen" => []) 
            for _ in 1:num_outages 
                component = rand(components)
                (component == "line") && (push!(scenarios[i]["branch"], rand(line_ids)))
                (component == "gen") && (push!(scenarios[i]["gen"], rand(gen_ids)))
            end 
            unique!(scenarios[i]["branch"])
            unique!(scenarios[i]["gen"])
        end 

        open(scenario_file, "w") do f 
            JSON.print(f, scenarios, 2) 
        end

        Base.run(`tar -zcvf $tar_file --absolute-paths $scenario_file`)
        Base.run(`rm -f $scenario_file`)
    else 
        df = DataFrame(CSV.File(cluster_file))
        gdf = groupby(df, :cluster_id)
        num_clusters = length(gdf)
        cluster_info = Dict()
        for i in 1:num_clusters
            cluster = gdf[i]
            buses = cluster[!, "bus_id"] |> collect 
            cluster_info[i] = buses
        end 
    end 
end 

main()