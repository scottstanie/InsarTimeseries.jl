using ArgParse

include("./InsarTimeseries.jl")
using .InsarTimeseries
include("./form_igrams.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--xlooks"
        default = 1
        arg_type = Int
        "--ylooks"
        default = 1
        arg_type = Int
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    collooks = parsed_args["xlooks"]
    rowlooks = parsed_args["ylooks"]
    @show rowlooks, collooks
    @time create_igrams(rowlooks, collooks)
end

main()
