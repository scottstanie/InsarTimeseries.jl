using ArgParse
include("./InsarTimeseries.jl")
using .InsarTimeseries

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--stack-average", "-s"
            help = "(Average a subset of independent igrams to get a linear solution."*
                   " If not, SBAS is used." 
            action = :store_true
        "--constant-velocity", "-c"
            help = "(SBAS) Use constant velocity (linear solution) to invert"
            action = :store_true
        "--outfile", "-o"
            default = "deformation.h5"
            help = "Name output .h5 file to save solution stack"
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
            # range_tester = x-> (x>0)
            help = "(SBAS) Strength of Tikhonov regularization"
        "--ref-station"
            arg_type = String
            help = "Name of a gps reference station to shift the stack to."*
                   " Will overwrite the current shifted dataset if specified"
        "path"
            help = "Path to directory with .h5 files"
            default = "."
            required = false
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
    ref_station = parsed_args["ref-station"]
    if !isnothing(ref_station)
        @time InsarTimeseries.shift_unw_file(unw_file, ref_station=ref_station, overwrite=true)
    end

    InsarTimeseries.run_inversion(unw_file,
                                  outfile=parsed_args["outfile"],
                                  use_stackavg=parsed_args["stack-average"],
                                  constant_velocity=parsed_args["constant-velocity"],
                                  ignore_geo_file=parsed_args["ignore-geo-file"],
                                  alpha=parsed_args["alpha"])
end


