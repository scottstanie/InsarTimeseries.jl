include("../src/InsarTimeseries.jl")
using .InsarTimeseries
using Test

@time begin
    @testset "Common" begin
        include("test_common.jl")
    end
    @testset "Stack averaging" begin
        include("test_stackavg.jl")
    end
    @testset "SBAS inversion" begin
        include("test_sbas.jl")
    end
end
