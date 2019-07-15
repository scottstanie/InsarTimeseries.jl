# TODO: is there a way to have build_executable import the cli.jl, but also 
# have the `main()` at the bottom? it errors now when it tries to run it....
# This feels kinds dumb now
include("./cli.jl")

main()
