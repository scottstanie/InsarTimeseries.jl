include("../src/InsarTimeseries.jl")
using .InsarTimeseries
using Test

@time begin
    @testset "Stack averaging" begin include("test_stackavg.jl") end
    @testset "SBAS inversion" begin include("test_sbas.jl") end
    @testset "Loading" begin include("test_loading.jl") end
end

