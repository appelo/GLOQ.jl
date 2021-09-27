#using GLOQ
using Test
using Random
using Plots
using DifferentialEquations
using LinearAlgebra
pyplot()
include("../src/GLOQ.jl")

@testset "GLOQ.jl" begin
    GLOQ.hello_world()
    include("forward_solve_constant_lindblad.jl")
    # Write your tests here.
end
