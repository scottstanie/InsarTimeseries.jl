using ArgParse

include("./InsarTimeseries.jl")
using .InsarTimeseries


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--overwrite"
            action = :store_true
            help = "Erase previous processed stack files"
        "--unw-stack-file"
            default = "unw_stack.h5"
            help = "name of .h5 file with unwrapped data"
        "--mask-stack-file"
            default = "masks.h5"
            help = "name of .h5 file with mask data"
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
            default = 5
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
    stack_flat_dset = STACK_FLAT_DSET
    stack_flat_shifted_dset = STACK_FLAT_SHIFTED_DSET

    @time InsarTimeseries.deramp_unw_file(unw_file, overwrite=false, stack_flat_dset=stack_flat_dset)

    ref_station = parsed_args["ref-station"]
    ref_row = parsed_args["ref-row"]
    ref_col = parsed_args["ref-col"]
    if !isnothing(ref_station) || !(isnothing(ref_row) || isnothing(ref_col))
        window = parsed_args["window"]
        @time InsarTimeseries.shift_unw_file(unw_file, ref_station=ref_station,
                                             ref_row=ref_col, ref_col=ref_col, overwrite=true, 
                                             window=window, stack_flat_dset=stack_flat_dset)
    end

end


