using ArgParse

include("./InsarTimeseries.jl")
using .InsarTimeseries


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--overwrite"
            action = :store_true
            help = "Erase previous processed stack files"
        "--igram-path"
            default = "."
            help = "Path containing .int files"
        "--geo-path"
            default = "../"
            help = "Path containing .geo files"
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
            range_tester = x-> (x>=1)
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

    overwrite = parsed_args["overwrite"]
    igram_path = parsed_args["igram-path"]
    geo_path = parsed_args["geo-path"]
    ref_station = parsed_args["ref-station"]
    ref_row = parsed_args["ref-row"]
    ref_col = parsed_args["ref-col"]
    window = parsed_args["window"]

    @time InsarTimeseries.prepare_stacks(igram_path, overwrite=overwrite, geo_path=geo_path, 
                                         ref_row=ref_row, ref_col=ref_col, ref_station=ref_station, 
                                         window=window)

end


main()
