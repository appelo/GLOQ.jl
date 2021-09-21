#using GLOQ
using Test
include("../src/GLOQ.jl")

@testset "GLOQ.jl" begin
    GLOQ.hello_world()
    include("forward_solve_constant_lindblad.jl")
    # Write your tests here.
end
