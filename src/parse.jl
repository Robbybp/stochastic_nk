using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--file", "-f"
        help = "pglib case file name with additional fields for n-k"
        arg_type = String
        default = "pglib_opf_case14_ieee_nk.m"
        
        "--path", "-p"
        help = "data directory path"
        arg_type = String
        default = "../data/"

        "--algo", "-a"
        help = "algorithm to run (full/Lshaped/Lshapedreg)"
        arg_type = String
        default = "Lshaped"

        "--timeout", "-t"
        help = "time limit for the run in seconds"
        arg_type = Int
        default = 3600

        "--gap", "-g"
        help = "optimality gap in %"
        arg_type = Float64
        default = 1e-3

        "--batchsize", "-s"
        help = "batch size for Lshaped algorithm"
        arg_type = Int
        default = 50

        "--batchid", "-b"
        help = "batch id for running the (id) batch Lshaped algorithm"
        arg_type = Int
        default = 1

        "--numbatches", "-n"
        help = "total number of batches in the SAA approximation"
        arg_type = Int
        default = 1

        "--interdict", "-i"
        help = "interdict lines or (lines and generators) (l/lg)"
        arg_type = String
        default = "lg"
    end

    return parse_args(s)
end
