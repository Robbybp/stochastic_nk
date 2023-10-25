""" data class that holds a solution to the problem """ 
struct SolutionDeterministic 
    lines::Vector{Int}
    generators::Vector{Int}
    load_shed::Float64 
    stats::Dict{Symbol,Any}
end 

SolutionDeterministic() = SolutionDeterministic([], [], NaN, Dict())

function Base.show(io::IO, solution::SolutionDeterministic)
    longest_field_name =
        maximum([
            length(string(fname)) for fname in fieldnames(SolutionDeterministic)
        ]) + 2
    printstyled(io, "**************** Solution ****************\n", color=:cyan)
    for name in fieldnames(SolutionDeterministic)
        sname = string(name)
        pname = sname * repeat(" ", longest_field_name - length(sname))
        if getfield(solution, name) === nothing
            println(io, pname, ": NA")
        else
            println(io, pname, ": ", getfield(solution, name))
        end
    end
    printstyled(io, "*****************************************\n", color=:cyan)
end 

""" data class that holds the results for the deterministic problem """ 
struct ResultsDeterministic
    num_iterations::Int 
    objective_value::Float64
    bound::Float64 
    run_time_in_seconds::Float64 
    optimality_gap::Float64 
    solution::SolutionDeterministic
end 

ResultsDeterministic() = ResultsDeterministic(0, NaN, NaN, NaN, NaN, SolutionDeterministic())
