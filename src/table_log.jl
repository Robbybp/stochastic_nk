function get_table_config(;problem=:deterministic)
    if problem == :stochastic
        # TODO: decide what to log for the stochastic N-k problem and fill in the fields
        fields = []
        field_chars = []
    else 
        fields = [:Iter, :LB, :UB, :Gap, :Time]
        field_chars = [10, 10, 10, 20, 10]
    end 
    return fields, field_chars
end

function print_table_header(fields, field_chars)
    println(repeat("=", sum(field_chars)))
    println(get_table_header_line(fields, field_chars))
    return println(repeat("=", sum(field_chars)))
end

function print_table_footer(fields, field_chars; symbol = "=")
    return println(repeat(symbol, sum(field_chars)))
end 

function get_table_header_line(fields, field_chars)
    ln = ""
    i = 1
    for f in fields
        f = string(f)
        padding = field_chars[i] - length(f)
        ln *= repeat(" ", trunc(Int, floor(padding / 2)))
        ln *= f
        ln *= repeat(" ", trunc(Int, ceil(padding / 2)))
        i += 1
    end
    return ln
end

function is_table_diff(fields, last_arr, new_arr)
    if length(last_arr) != length(new_arr)
        return true
    end

    time_idx = findfirst(fields .== :Time)
    if time_idx !== nothing
        last_arr = vcat(last_arr[1:time_idx-1], last_arr[time_idx+1:end])
        new_arr = vcat(new_arr[1:time_idx-1], new_arr[time_idx+1:end])
    end
    for i in eachindex(last_arr)
        last_arr[i] != new_arr[i] && return true
    end
    return false
end

function print_table_cp(iter, lb, ub, gap, start_time, fields, field_chars; last_arr = [])
    table_line, table_arr = get_table_line_cp(iter, lb, ub, gap, start_time, fields, field_chars)
    is_table_diff(fields, last_arr, table_arr) && println(table_line)
    return table_arr
end

function get_table_line_cp(iter, lb, ub, gap, start_time, fields, field_chars)
    arr = []
    ln = ""
    precision = 4

    for (i, f) in enumerate(fields)
        val = ""
        if f == :Iter
            val = string(Int(iter))
        elseif f == :LB
            val = isnan(lb) ? "-" : string(round(lb; digits = precision))
        elseif f == :UB
            val = isnan(ub) ? "-" : string(round(ub; digits = precision))
        elseif f == :Time
            val = string(Int(ceil((now() - start_time).value/1000.0))) * "s"
        elseif f == :Gap 
            val = isinf(gap) ? "Inf" :  "⫹" * string(round(gap; digits = 2)) * "%"
        end
        if length(val) > field_chars[i]
            # too long to display shouldn't happen normally but is better than error
            # if it happens
            val = "t.l."
        end

        padding = field_chars[i] - length(val)
        ln *= repeat(" ", trunc(Int, floor(padding / 2)))
        ln *= val
        ln *= repeat(" ", trunc(Int, ceil(padding / 2)))
        push!(arr, val)
    end
    return ln, arr
end