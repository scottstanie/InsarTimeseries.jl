using ArgParse

include("./InsarTimeseries.jl")
using .InsarTimeseries

# # TODO: make this work for single arg
# import ArgParse.parse_item
# function ArgParse.parse_item(it_type::Tuple{Int, Int}, x::AbstractString) 
#     tup = split(replace(x, " "=>""), ",")
#     return convert(it_type, tup)
# end

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--stack-average", "-s"
            action = :store_true
            help = "(Average a subset of independent igrams to get a linear solution."*
                   " If not, SBAS is used." 
        "--constant-velocity", "-c"
            action = :store_true
            help = "(SBAS) Use constant velocity (linear solution) to invert"
        "--outfile", "-o"
            default = "deformation.h5"
            help = "Name output .h5 file to save solution stack"
        "--max-temporal-baseline"
            arg_type = Int
            range_tester = x-> (x>0)
            help = "Max span allowed between igrams"
        "--ignore-geo-file"
            help = "name of file with list of .geo  dates to ignore"
        "--unw-stack-file"
            default = "unw_stack.h5"
            help = "name of .h5 file with unwrapped data"
        "--mask-stack-file"
            default = "masks.h5"
            help = "name of .h5 file with mask data"
        "--alpha"
            default = 0.0f0
            arg_type = Float32
            range_tester = x-> (x>=0.0f0)
            help = "(SBAS) Strength of Tikhonov regularization"
        "--L1"
            action = :store_true
            default = false
            help = "Use L1 norm for SBAS cost instead of L2 least squares"
        # TODO: do i care to add this path option? 
        # "path"
        #     help = "Path to directory with .h5 files"
        #     default = "."
        #     required = false
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        @show arg, val
    end
    unw_file = parsed_args["unw-stack-file"]

    # Assign the appropriate dataset name based on the order
    InsarTimeseries.run_inversion(unw_file,
                                  outfile=parsed_args["outfile"],
                                  use_stackavg=parsed_args["stack-average"],
                                  constant_velocity=parsed_args["constant-velocity"],
                                  ignore_geo_file=parsed_args["ignore-geo-file"],
                                  max_temporal_baseline=parsed_args["max-temporal-baseline"],
                                  alpha=parsed_args["alpha"],
                                  L1=parsed_args["L1"])
end


