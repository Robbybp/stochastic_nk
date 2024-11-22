""" data class that holds a solution to the problem """ 
struct Solution
    lines::Vector{Int}
    generators::Vector{Int}
    buses::Vector{Int}
    load_shed::Float64 
    stats::Dict{Symbol,Any}
end 

Solution() = Solution([], [], [], NaN, Dict())

function Base.show(io::IO, solution::Solution)
    longest_field_name =
        maximum([
            length(string(fname)) for fname in fieldnames(Solution)
        ]) + 2
    printstyled(io, "\n**************** Solution ****************\n", color=:cyan)
    for name in fieldnames(Solution)
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

""" data class that holds the results """ 
struct Results
    num_iterations::Int 
    objective_value::Float64
    bound::Float64 
    run_time_in_seconds::Float64 
    optimality_gap::Float64 
    solution::Solution
end 

Results() = Results(0, NaN, NaN, NaN, NaN, Solution())
