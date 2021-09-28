using Random
using Plots
using DifferentialEquations,DiffEqFlux
using Optim
using ForwardDiff
include("../src/GLOQ.jl")
pyplot()
#import GLOQ
############################################################
# Step 1: Generate Synthetic Data
############################################################
HK = [0.0 0.0 1.0;
      0.0 0.0 0.0;
      1.0 0.0 0.0]
HS = [0.0 1.0 0.0;
     -1.0 0.0 0.0;
      0.0 0.0 0.0]
L1 = [0.0 0.25 0.0;
      0.0 0.0  0.5;
      0.0 0.0  0.0]
L2 = [0.0 0.0 0.0;
      0.0 0.5 0.0;
      0.0 0.0 1.0]
L_list = (L1,L2)
LK,LS,LD = GLOQ.make_lindblad_operator(HK,HS,L_list)
println("Lindblad operator assembled")
# define the time interval
tspan_length = 21
t_span = range(0.0,pi,length=tspan_length)
# forward solves with different initial conditions,
# use exponential_solver
u1 = [1.0;0;0]
rho1_real,rho1_imag_adj = GLOQ.exponential_solver(real(u1),-imag(u1),LK,LS,LD,t_span;initial_type = "states")
u2 = [0;1.0;0]
rho2_real,rho2_imag_adj = GLOQ.exponential_solver(real(u2),-imag(u2),LK,LS,LD,t_span;initial_type = "states")
u3 = [1.0/sqrt(2.0);-im/sqrt(2.0);0]
rho3_real,rho3_imag_adj = GLOQ.exponential_solver(real(u3),-imag(u3),LK,LS,LD,t_span;initial_type = "states")
u4 = [0.0;1.0/sqrt(2.0);-im/sqrt(2.0)]
rho4_real,rho4_imag_adj = GLOQ.exponential_solver(real(u4),-imag(u4),LK,LS,LD,t_span;initial_type = "states")
u5 = [1.0/sqrt(2.0);0.0;-im/sqrt(2.0)]
rho5_real,rho5_imag_adj = GLOQ.exponential_solver(real(u5),-imag(u5),LK,LS,LD,t_span;initial_type = "states")

population1_history = GLOQ.get_population(rho1_real)
population2_history = GLOQ.get_population(rho2_real)
population3_history = GLOQ.get_population(rho3_real)
population4_history = GLOQ.get_population(rho4_real)
population5_history = GLOQ.get_population(rho5_real)

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
rho1_initial_real,rho1_initial_adj_imag = GLOQ.convert_state_to_density(real(u1),-imag(u1))
rho2_initial_real,rho2_initial_adj_imag = GLOQ.convert_state_to_density(real(u2),-imag(u2))
rho3_initial_real,rho3_initial_adj_imag = GLOQ.convert_state_to_density(real(u3),-imag(u3))
rho4_initial_real,rho4_initial_adj_imag = GLOQ.convert_state_to_density(real(u4),-imag(u4))
rho5_initial_real,rho5_initial_adj_imag = GLOQ.convert_state_to_density(real(u5),-imag(u5))
rho1_initial = [rho1_initial_real;rho1_initial_adj_imag]
rho2_initial = [rho2_initial_real;rho2_initial_adj_imag]
rho3_initial = [rho3_initial_real;rho3_initial_adj_imag]
rho4_initial = [rho4_initial_real;rho4_initial_adj_imag]
rho5_initial = [rho5_initial_real;rho5_initial_adj_imag]
# p[1]-p[5] upper diagonal of the Hamiltonian part
# p[6]-p[9] parameters in the Lindblad term
function loss(p)
      _HK = [0  0  p[2];
             0  0  0;
             p[2] 0  0]
      _HS = [0  p[1]  0;
            -p[1]  0  0;
             0  0  0]
      _L1 = [0.0 p[3] 0.0;
            0.0 0.0  p[4];
            0.0 0.0  0.0]
      _L2 = [0.0 0.0 0.0;
            0.0 p[5] 0.0;
            0.0 0.0 2.0*p[5]]
      _L_list = (_L1,_L2)
      _LK,_LS,_LD = GLOQ.make_lindblad_operator(_HK,_HS,_L_list)

      _problem = GLOQ.LindbladODEProblem(rho1_initial_real,rho1_initial_adj_imag,
                                        _LK,_LS,_LD,
                                         t_span[end];initial_type="density")
      rho1_data = solve(_problem,saveat=t_span)
      rho2_data = solve(remake(_problem, u0=rho2_initial), saveat=t_span)
      rho3_data = solve(remake(_problem, u0=rho3_initial), saveat=t_span)
      rho4_data = solve(remake(_problem, u0=rho4_initial), saveat=t_span)
      rho5_data = solve(remake(_problem, u0=rho5_initial), saveat=t_span)

      pop_data1 = GLOQ.get_population(Array(rho1_data[1:Int64(end/2),:]))
      pop_data2 = GLOQ.get_population(Array(rho2_data[1:Int64(end/2),:]))
      pop_data3 = GLOQ.get_population(Array(rho3_data[1:Int64(end/2),:]))
      pop_data4 = GLOQ.get_population(Array(rho4_data[1:Int64(end/2),:]))
      pop_data5 = GLOQ.get_population(Array(rho5_data[1:Int64(end/2),:]))

      _loss = sum(abs2, pop_data1-noisy_data1)
      _loss += sum(abs2, pop_data2-noisy_data2)
      _loss += sum(abs2, pop_data3-noisy_data3)
      _loss += sum(abs2, pop_data4-noisy_data4)
      _loss += sum(abs2, pop_data5-noisy_data5)
      return _loss
end

p_initial = [1.0;1.0;0.5;0.5;0.5]
loss_initial = loss(p_initial)
println("Initial loss: ",loss_initial)
p_true = [1.0 1.0 0.25 0.5 0.5]

@time p_DEFlux = DiffEqFlux.sciml_train(loss, p_initial)
@time p_Optim = Optim.optimize(loss, p_initial)
println("True p: ",p_true)
println("DiffEqFlux p: ",p_DEFlux.u)
println("DiffEqFlux training history: ",p_DEFlux)
println("Optim p: ",p_Optim.minimizer)
println("Optim training history",p_Optim)


################################
function loss3(p)
      _HK = [0    0    p[2];
             0    p[3] p[4];
             p[2] p[4] p[5]]
      _HS = [0  p[1]  0;
            -p[1]  0  0;
             0  0  0]
      _L1 = [0.0 p[6] 0.0;
            0.0 0.0  p[7];
            0.0 0.0  0.0]
      _L2 = [0.0 0.0 0.0;
            0.0 p[8] 0.0;
            0.0 0.0 2.0*p[8]]
      _L_list = (_L1,_L2)
      _LK,_LS,_LD = GLOQ.make_lindblad_operator(_HK,_HS,_L_list)

      _problem = GLOQ.LindbladODEProblem(rho1_initial_real,rho1_initial_adj_imag,
                                              _LK,_LS,_LD,
                                               t_span[end];initial_type="density")
      rho1_data = solve(_problem,saveat=t_span)
      rho2_data = solve(remake(_problem, u0=rho2_initial), saveat=t_span)
      rho3_data = solve(remake(_problem, u0=rho3_initial), saveat=t_span)
      rho4_data = solve(remake(_problem, u0=rho4_initial), saveat=t_span)
      rho5_data = solve(remake(_problem, u0=rho5_initial), saveat=t_span)

      pop_data1 = GLOQ.get_population(Array(rho1_data[1:Int64(end/2),:]))
      pop_data2 = GLOQ.get_population(Array(rho2_data[1:Int64(end/2),:]))
      pop_data3 = GLOQ.get_population(Array(rho3_data[1:Int64(end/2),:]))
      pop_data4 = GLOQ.get_population(Array(rho4_data[1:Int64(end/2),:]))
      pop_data5 = GLOQ.get_population(Array(rho5_data[1:Int64(end/2),:]))

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

##########################################
function test_auto_diff(p)
      _HK = [0  0  p[2];
             0  0  0;
             p[2] 0  0]
      _HS = [0  p[1]  0;
            -p[1]  0  0;
             0  0  0]
      _L1 = [0.0 p[3] 0.0;
            0.0 0.0  p[4];
            0.0 0.0  0.0]
      _L2 = [0.0 0.0 0.0;
            0.0 p[5] 0.0;
            0.0 0.0 2.0*p[5]]
      _L_list = (_L1,_L2)
      _LK,_LS,_LD = GLOQ.make_lindblad_operator(_HK,_HS,_L_list)
      _problem = GLOQ.LindbladODEProblem(rho1_initial_real,rho1_initial_adj_imag,
                                        _LK,_LS,_LD,
                                         t_span[end];initial_type="density")
      rho1_data = solve(_problem)
#      return norm(Array(rho1_data[:,end]),2)^2 # AD works
      pop_data1 = transpose(rho1_data[1:4:end,:])
      return norm(pop_data1[:,end],2)
end
@time Zygote.gradient(test_auto_diff, p_true)
println("auto-diff done")
