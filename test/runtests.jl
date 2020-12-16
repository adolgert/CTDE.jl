using CTDE
using Test
using SafeTestsets

@testset "CTDE.jl" begin
    include("test_prefixsearch.jl")
    include("test_direct.jl")
end
