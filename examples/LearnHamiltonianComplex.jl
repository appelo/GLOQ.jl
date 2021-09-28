using Random
using Plots
using DifferentialEquations,DiffEqFlux
using Optim
include("../src/GLOQ.jl")
pyplot()
#import GLOQ
############################################################
# Step 1: Generate Synthetic Data
############################################################
H = [ 0 im 1.0;
     -im 0 0;
     1.0   0 0]
L1 = [0.0 0.25 0.0;
      0.0 0.0  0.5;
      0.0 0.0  0.0]
L2 = [0.0 0.0 0.0;
      0.0 0.5 0.0;
      0.0 0.0 1.0]
L_list = (L1,L2)
LH,LD = GLOQ.make_lindblad_operator_complex(H,L_list)
L = LH+LD
println("Lindblad operator assembled")
# define the time interval
tspan_length = 21
t_span = range(0.0,pi,length=tspan_length)
# forward solves with different initial conditions,
# use exponential_solver
u1 = [1.0;0;0]
rho1_history = GLOQ.exponential_solver_complex(u1,L,t_span;initial_type = "states")
u2 = [0;1.0;0]
rho2 = (u2*u2')[:]
rho2_history = GLOQ.exponential_solver_complex(u2,L,t_span;initial_type = "states")
u3 = [1.0/sqrt(2.0);-im/sqrt(2.0);0]
rho3_history = GLOQ.exponential_solver_complex(u3,L,t_span;initial_type = "states")
u4 = [0.0;1.0/sqrt(2.0);-im/sqrt(2.0)]
rho4_history = GLOQ.exponential_solver_complex(u4,L,t_span;initial_type = "states")
u5 = [1.0/sqrt(2.0);0.0;-im/sqrt(2.0)]
rho5_history = GLOQ.exponential_solver_complex(u5,L,t_span;initial_type = "states")

population1_history = GLOQ.get_population(rho1_history)
population2_history = GLOQ.get_population(rho2_history)
population3_history = GLOQ.get_population(rho3_history)
population4_history = GLOQ.get_population(rho4_history)
population5_history = GLOQ.get_population(rho5_history)

fig = plot(t_span,population1_history,lab=:false)
plot!(fig,t_span,population2_history,lab=:false)
plot!(fig,t_span,population3_history,lab=:false)
plot!(fig,t_span,population4_history,lab=:false)
plot!(fig,t_span,population5_history,lab=:false)
display(fig)

# add noise to the synthetic data
noise_amp = 0.025
noisy_data1 = population1_history+noise_amp.*rand(Float64,size(population1_history))
noisy_data2 = population2_history+noise_amp.*rand(Float64,size(population2_history))
noisy_data3 = population3_history+noise_amp.*rand(Float64,size(population3_history))
noisy_data4 = population4_history+noise_amp.*rand(Float64,size(population4_history))
noisy_data5 = population5_history+noise_amp.*rand(Float64,size(population5_history))
# present the noisy data
fig = plot(t_span,noisy_data1,lab=:false)
plot!(fig,t_span,noisy_data2,lab=:false)
plot!(fig,t_span,noisy_data3,lab=:false)
plot!(fig,t_span,noisy_data4,lab=:false)
plot!(fig,t_span,noisy_data5,lab=:false)
display(fig)
############################################################
# Step 2: Learn the operators
############################################################
rho1 = convert(Vector{ComplexF64},(u1*u1')[:])
rho2 = convert(Vector{ComplexF64},(u2*u2')[:])
rho3 = convert(Vector{ComplexF64},(u3*u3')[:])
rho4 = convert(Vector{ComplexF64},(u4*u4')[:])
rho5 = convert(Vector{ComplexF64},(u5*u5')[:])
# p[1]-p[5] upper diagonal of the Hamiltonian part
# p[6]-p[9] parameters in the Lindblad term
function loss(p)
      _X = [ 0 p[1]*im p[2];
             0   0  0;
             0   0  0]
      _H = 0.5.*(_X+adjoint(_X))
      _L1 = [0.0 p[3] 0.0;
            0.0 0.0  p[4];
            0.0 0.0  0.0]
      _L2 = [0.0 0.0 0.0;
            0.0 p[5] 0.0;
            0.0 0.0 2.0*p[5]]
      _L_list = (_L1,_L2)
      _LH,_LD = GLOQ.make_lindblad_operator_complex(_H,_L_list)
      _L = _LH+_LD

      _problem = GLOQ.LindbladODEProblemComplex(rho1,_L,t_span[end];initial_type="density")
      rho1_data = solve(_problem,saveat=t_span)
      rho2_data = solve(remake(_problem, u0=rho2), saveat=t_span)
      rho3_data = solve(remake(_problem, u0=rho3), saveat=t_span)
      rho4_data = solve(remake(_problem, u0=rho4), saveat=t_span)
      rho5_data = solve(remake(_problem, u0=rho5), saveat=t_span)

      pop_data1 = GLOQ.get_population(Array(rho1_data))
      pop_data2 = GLOQ.get_population(Array(rho2_data))
      pop_data3 = GLOQ.get_population(Array(rho3_data))
      pop_data4 = GLOQ.get_population(Array(rho4_data))
      pop_data5 = GLOQ.get_population(Array(rho5_data))

      _loss = sum(abs2, pop_data1-noisy_data1)
      _loss += sum(abs2, pop_data2-noisy_data2)
      _loss += sum(abs2, pop_data3-noisy_data3)
      _loss += sum(abs2, pop_data4-noisy_data4)
      _loss += sum(abs2, pop_data5-noisy_data5)
      return _loss
end

p_initial = [1.0;1.0;0.5;0.5;0.5]
#loss_initial = loss(p_initial)
#println(loss_initial)
p_true = [2.0 2.0 0.25 0.5 0.5]
@time p_DEFlux = DiffEqFlux.sciml_train(loss, p_initial)
@time p_Optim = Optim.optimize(loss, p_initial)
println("True p: ",p_true)
println("DEFlux p: ",p_DEFlux.u)
println("Optim p: ",p_Optim.minimizer)
println(p_DEFlux)
println(p_Optim)

################################
function loss3(p)
      _X = [ 0 p[1]*im p[2];
             0   p[3]  p[4];
             0   0     p[5]]
      _H = 0.5.*(_X+adjoint(_X))
      _L1 = [0.0 p[6] 0.0;
            0.0 0.0  p[7];
            0.0 0.0  0.0]
      _L2 = [0.0 0.0 0.0;
            0.0 p[8] 0.0;
            0.0 0.0 2.0*p[8]]
      _L_list = (_L1,_L2)
      _LH,_LD = GLOQ.make_lindblad_operator_complex(_H,_L_list)
      _L = _LH+_LD

      _problem = GLOQ.LindbladODEProblemComplex(rho1,_L,t_span[end];initial_type="density")
      rho1_data = solve(_problem,saveat=t_span)
      rho2_data = solve(remake(_problem, u0=rho2), saveat=t_span)
      rho3_data = solve(remake(_problem, u0=rho3), saveat=t_span)
      rho4_data = solve(remake(_problem, u0=rho4), saveat=t_span)
      rho5_data = solve(remake(_problem, u0=rho5), saveat=t_span)

      pop_data1 = GLOQ.get_population(Array(rho1_data))
      pop_data2 = GLOQ.get_population(Array(rho2_data))
      pop_data3 = GLOQ.get_population(Array(rho3_data))
      pop_data4 = GLOQ.get_population(Array(rho4_data))
      pop_data5 = GLOQ.get_population(Array(rho5_data))

      _loss = sum(abs2, pop_data1-noisy_data1)
      _loss += sum(abs2, pop_data2-noisy_data2)
      _loss += sum(abs2, pop_data3-noisy_data3)
      _loss += sum(abs2, pop_data4-noisy_data4)
      _loss += sum(abs2, pop_data5-noisy_data5)
      return _loss
end

p3_initial = [1.0;1.0;0.25;0.1;0.1;0.5;0.5;0.5]
@time p3_DEFlux = DiffEqFlux.sciml_train(loss3, p3_initial)
@time p3_Optim = Optim.optimize(loss3, p3_initial)
println("True p: ",p_true)
println("DEFlux p: ",p3_DEFlux.u)
println("Optim p: ",p3_Optim.minimizer)
println(p3_DEFlux)
println(p3_Optim)
################################################
#=
function loss2(p)
      _H = p[1]+1.0im.*p[2]
      _L1 = p[3]
      _L2 = p[4]
      _L_list = (_L1,_L2)
      _LH,_LD = GLOQ.make_lindblad_operator(_H,_L_list)
      _L = _LH+_LD

      _problem = GLOQ.LindbladODEProblem(rho1,_L,t_span[end];initial_type="density")
      rho1_data = solve(_problem,saveat=t_span)
      rho2_data = solve(remake(_problem, u0=rho2), saveat=t_span)
      rho3_data = solve(remake(_problem, u0=rho3), saveat=t_span)
      rho4_data = solve(remake(_problem, u0=rho4), saveat=t_span)
      rho5_data = solve(remake(_problem, u0=rho5), saveat=t_span)

      pop_data1 = GLOQ.get_population(Array(rho1_data))
      pop_data2 = GLOQ.get_population(Array(rho2_data))
      pop_data3 = GLOQ.get_population(Array(rho3_data))
      pop_data4 = GLOQ.get_population(Array(rho4_data))
      pop_data5 = GLOQ.get_population(Array(rho5_data))

      _loss = sum(abs2, pop_data1-noisy_data1)
      _loss += sum(abs2, pop_data2-noisy_data2)
      _loss += sum(abs2, pop_data3-noisy_data3)
      _loss += sum(abs2, pop_data4-noisy_data4)
      _loss += sum(abs2, pop_data5-noisy_data5)
      return _loss
end

Hr_initial = [ 0 0 1.2;
               0 0 0;
            1.2  0 0]
Hi_initial = [ 0 1.5 0;
              -1.5 0 0;
               1  0 0]
L1_initial = [0.0 0.5 0.0;
              0.0 0.0 0.3;
              0.0 0.0 0.0]
L2_initial = [0.0 0.0 0.0;
              0.0 0.3 0.0;
              0.0 0.0 0.6]
p2_initial = (Hr_initial,Hi_initial,L1_initial,L2_initial)
loss2_initial = loss2(p2_initial)

@time p_opt2 = DiffEqFlux.sciml_train(loss, p2_initial,BFGS(),maxiters=100)
#println("True p: ",p_true)
println("Optimized Operators: ",p_opt2.u)
println(H," ",L1," ",L2)
=#
