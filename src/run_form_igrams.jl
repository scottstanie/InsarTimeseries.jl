using ArgParse

include("./InsarTimeseries.jl")
using .InsarTimeseries
include("./form_igrams.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--looks"
            default = 1
            arg_type = Int
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    looks = parsed_args["looks"]
    @show looks
    @time create_igrams(looks, looks)

main()
