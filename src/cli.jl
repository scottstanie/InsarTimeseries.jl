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
            help = "(Average a subset of independent igrams to get a linear solution."*
                   " If not, SBAS is used." 
            action = :store_true
        "--constant-velocity", "-c"
            help = "(SBAS) Use constant velocity (linear solution) to invert"
            action = :store_true
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
        "--ref-station"
            arg_type = String
            help = "Name of a gps reference station to shift the stack to."*
                   " Will overwrite the current shifted dataset if specified"
        "--ref-row"
            arg_type = Int
        "--ref-col"
            arg_type = Int
        "--window"
            arg_type = Int
            range_tester = x-> (x>0)
            help = "Size of window to use to shift stack"
            default = 3
        "--order"
            arg_type = Int
            range_tester = x -> (x == 1 || x == 2)
            default = 1
            help = "Order of ramp to use for deramping .unw files"
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

    order = parsed_args["order"]

    # Assign the appropriate dataset name based on the order
    stack_flat_dset = (order == 1) ? InsarTimeseries.STACK_FLAT_DSET1 : InsarTimeseries.STACK_FLAT_DSET2
    stack_flat_shifted_dset = (order == 1) ? InsarTimeseries.STACK_FLAT_SHIFTED_DSET1 : InsarTimeseries.STACK_FLAT_SHIFTED_DSET2

    @time InsarTimeseries.deramp_unw_file(unw_file, order=order, overwrite=false, stack_flat_dset=stack_flat_dset)

    ref_station = parsed_args["ref-station"]
    ref_row = parsed_args["ref-row"]
    ref_col = parsed_args["ref-col"]
    if !isnothing(ref_station) || !(isnothing(ref_row) || isnothing(ref_col))
        window = parsed_args["window"]
        @time InsarTimeseries.shift_unw_file(unw_file, ref_station=ref_station,
                                             ref_row=ref_col, ref_col=ref_col, overwrite=true, 
                                             window=window, stack_flat_dset=stack_flat_dset)
    end

    InsarTimeseries.run_inversion(unw_file,
                                  outfile=parsed_args["outfile"],
                                  use_stackavg=parsed_args["stack-average"],
                                  constant_velocity=parsed_args["constant-velocity"],
                                  ignore_geo_file=parsed_args["ignore-geo-file"],
                                  max_temporal_baseline=parsed_args["max-temporal-baseline"],
                                  alpha=parsed_args["alpha"])
end


