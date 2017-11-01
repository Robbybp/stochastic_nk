
function create_master_model(scenarios, ref::Dict{Symbol,Any}, config::Dict{String,Any}, model=Model())
    
    ref = ref[:nw][0]
    numscenarios = config["batchsize"]
    config["budget"] == sum(scenarios[1,:])

    # interdiction variables
    @variable(model, x[i in keys(ref[:branch])], Bin)
    @variable(model, y[i in keys(ref[:gen])], Bin)

    # lifted variables for multi-cut Lshaped
    @variable(model, θ[1:numscenarios] >= 0)

    # master constraints
    @constraint(model, sum(x) + sum(y) == config["budget"])
    (config["interdict"] == "l") && @constraint(model, sum(y) == 0)

    return model

end
