# TODO: figure out if there's a way to get around this hacky file-splitting thing
include("./cli.jl")


Base.@ccallable function julia_main(ARGS::Vector{String})::Cint
    main()
    return 0
end

