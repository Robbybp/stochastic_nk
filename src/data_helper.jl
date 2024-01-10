""" set generation limits to zero when pmax/qmax are positive """ 
function modify_generation_limits(case_data::Dict)
    for (_, gen) in get(case_data, "gen", [])
        (gen["pmax"] > 0.0) && (gen["pmin"] = 0.0)
        (gen["qmax"] > 0.0) && (gen["qmin"] = 0.0)
    end 
end

""" modify phase angle bounds """ 
function modify_phase_angle_bounds(case_data::Dict)
    for (_, branch) in get(case_data, "branch", [])
        branch["angmin"] = -pi 
        branch["angmax"] = pi
    end 
end 

""" add total load to case_data """ 
function add_total_load_info(case_data::Dict)
    pd = sum([load["pd"] for (_, load) in case_data["load"]]; init=0.0)
    ps = sum([shunt["gs"] for (_, shunt) in case_data["shunt"]]; init=0.0)
    case_data["total_load"] = pd + ps
end 