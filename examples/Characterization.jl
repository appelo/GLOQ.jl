using Zygote
using Random
using DelimitedFiles
using LinearAlgebra
using GalacticOptim,Optim,NLopt#,BlackBoxOptim
using DifferentialEquations#,DiffEqFlux
#using ForwardDiff
using Plots
using LaTeXStrings
using GLOQ
#using ReverseDiff
#include("../src/GLOQ.jl")
pyplot()
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

N_states = 4;
freqs = [4.10817777; 3.88303940; 3.61937188]
omegas = 2.0*pi.*freqs
gamma1   = [2.59451982e-05; 4.54296537e-05; 0.0]
gamma2   = [2.58005604e-05; 9.00000000e-05; 0.0]

# driving frequency
omr_t1_01 = 2.0*pi*4.1081786
omr_t1_12 = 2.0*pi*3.8830355
omr_ramsey_01 = 2.0*pi*(4.1081786-5e-4)
omr_ramsey_12 = 2.0*pi*(3.8830355-2.0e-3)
omr_echo_01 = 2.0*pi*4.1081786
omr_echo_12 = 2.0*pi*3.8830355

# charge noise
char_noise = 2.0*pi.*[-1.017e-5,1.058e-4,0]

# width of the control signal
width = 21.25
TC_01 = 2.5*width

width = 17.00
TC_12 = 2.5*width

# in nano second
T_T1_01 = 200.0*1000.0
T_T1_12 = 100.0*1000.0
T_Ramsey_01 = 40.0*1000.0
T_Ramsey_12 = 20.0*1000.0
T_Echo_01 = 40.0*1000.0
T_Echo_12 = 20.0*1000.0

# data size used for characterization
n1use = 201
n2use = 201
n3use = 201
n4use = 281
n5use = 201
n6use = 201

# time points evaluated
dt_t1_01 = T_T1_01/200
dt_t1_12 = T_T1_12/200
dt_ramsey_01 = T_Ramsey_01/200
dt_ramsey_12 = T_Ramsey_12/400
dt_echo_01 = T_Echo_01/200
dt_echo_12 = T_Echo_12/200

t_t1_01 = zeros(Float64,n1use)
for i = 1:n1use
	t_t1_01[i] = (i-1)*dt_t1_01
end
T_evaluated_t1_01 = t_t1_01[end]

t_t1_12 = zeros(Float64,n2use)
for i = 1:n2use
	t_t1_12[i] = (i-1)*dt_t1_12
end
T_evaluated_t1_12 = t_t1_12[end]

t_ramsey_01 = zeros(Float64,n3use)
for i = 1:n3use
	t_ramsey_01[i] = (i-1)*dt_ramsey_01
end
T_evaluated_ramsey_01 = t_ramsey_01[end]

t_ramsey_12 = zeros(Float64,n4use)
for i = 1:n4use
	t_ramsey_12[i] = (i-1)*dt_ramsey_12
end
T_evaluated_ramsey_12 = t_ramsey_12[end]

t_echo_01 = zeros(Float64,n5use)
for i = 1:n5use
	t_echo_01[i] = (i-1)*dt_echo_01
end
T_evaluated_echo_01 = t_echo_01[end]

t_echo_12 = zeros(Float64,n6use)
for i = 1:n6use
	t_echo_12[i] = (i-1)*dt_echo_12
end
T_evaluated_echo_12 = t_echo_12[end]

data_t1_01 = readdlm("data/CD8_Alice_0-1_T1.txt",)
data_t1_01 = transpose(data_t1_01[:,1:n1use])
data_t1_12 = readdlm("data/CD8_Alice_1-2_T1.txt")
data_t1_12 = transpose(data_t1_12[:,1:n2use])
data_ramsey_01 = readdlm("data/Alice_ramsey_0-1_20200722.txt")
data_ramsey_01 = transpose(data_ramsey_01[:,1:n3use])
data_ramsey_12 = readdlm("data/Alice_ramsey_1-2_20200723.txt")
data_ramsey_12 = transpose(data_ramsey_12[:,1:n4use])
data_echo_01 = readdlm("data/CD8_Alice_0-1_echo.txt")
data_echo_01 = transpose(data_echo_01[:,1:n5use])
data_echo_12 = readdlm("data/CD8_Alice_1-2_echo.txt")
data_echo_12 = transpose(data_echo_12[:,1:n6use])

#println(size(data_t1_01))
#println(size(data_t1_12))
#println(size(data_ramsey_01))
#println(size(data_ramsey_12))
#println(size(data_echo_01))
#println(size(data_echo_12))

# initial conditions
u0_t1_01 = [1.0;0.0;0.0;0.0]
u0_t1_12 = [0.0;1.0;0.0;0.0]
u0_ramsey_01 = [1.0;0.0;0.0;0.0]
u0_ramsey_12 = [0.0;1.0;0.0;0.0]
u0_echo_01 = [1.0;0.0;0.0;0.0]
u0_echo_12 = [0.0;1.0;0.0;0.0]

v0_t1_01 = [0.0;0.0;0.0;0.0]
v0_t1_12 = [0.0;0.0;0.0;0.0]
v0_ramsey_01 = [0.0;0.0;0.0;0.0]
v0_ramsey_12 = [0.0;0.0;0.0;0.0]
v0_echo_01 = [0.0;0.0;0.0;0.0]
v0_echo_12 = [0.0;0.0;0.0;0.0]

########################################################
# Define loss function
########################################################
function loss(p)#,p_keywords)
	_rho_Ramsey_u,_rho_Ramsey_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				     [p[1:2];omegas[3]],omr_ramsey_12,
					 [p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 initial_state,
					 TC_12,t_ramsey_12,N_states)
	_population_Ramsey12 = GLOQ.get_population(_rho_Ramsey_u)
	_loss = sum(abs2,_population_Ramsey12[:,1:3]-data_ramsey_12).*dt_ramsey_12/1000
	return _loss
end

global _function_call = 0
global _iteration_number = 0
function loss_gala(p,p_keywords)
	# Ramsey experiments
	_rho_Ramsey_01_u,_rho_Ramsey_01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
					 [2.0*pi.*p[1:2];omegas[3]],omr_ramsey_01,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 0,
					 TC_01,t_ramsey_01,N_states)
	_population_Ramsey_01 = GLOQ.get_population(_rho_Ramsey_01_u)

	_rho_Ramsey_12_u,_rho_Ramsey_12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				     [2.0*pi.*p[1:2];omegas[3]],omr_ramsey_12,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 1,
					 TC_12,t_ramsey_12,N_states)
	_population_Ramsey_12 = GLOQ.get_population(_rho_Ramsey_12_u)

	# Echo experiments
	_rho_Echo_01_u,_rho_Echo_01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
					 [2.0*pi.*p[1:2];omegas[3]],omr_echo_01,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 0,
					 TC_01,t_echo_01,N_states)
	_population_Echo_01 = GLOQ.get_population(_rho_Echo_01_u)

	_rho_Echo_12_u,_rho_Echo_12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
					 [2.0*pi.*p[1:2];omegas[3]],omr_echo_12,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 1,
					 TC_12,t_echo_12,N_states)
	_population_Echo_12 = GLOQ.get_population(_rho_Echo_12_u)

	# T1-decay experiments
	_rho_T1_01_u,_rho_T1_01_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
					 [2.0*pi.*p[1:2];omegas[3]],omr_t1_01,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 0,
					 TC_01,t_t1_01,N_states)
	_population_T1_01 = GLOQ.get_population(_rho_T1_01_u)

	_rho_T1_12_u,_rho_T1_12_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
					 [2.0*pi.*p[1:2];omegas[3]],omr_t1_12,
					 [2.0*pi.*p[7:8];char_noise[3]],
					 [p[3:4];0.0],[p[5:6];0.0],
					 1,
					 TC_12,t_t1_12,N_states)
	_population_T1_12 = GLOQ.get_population(_rho_T1_12_u)

	################################################################################
	_loss = sum(abs2,_population_Ramsey_01[:,1:3]-data_ramsey_01)*dt_ramsey_01/T_evaluated_ramsey_01+
			sum(abs2,_population_Ramsey_12[:,1:3]-data_ramsey_12)*dt_ramsey_12/T_evaluated_ramsey_12
			sum(abs2,_population_Echo_01[:,1:3]-data_echo_01).*dt_echo_01/T_evaluated_echo_01+
			sum(abs2,_population_Echo_12[:,1:3]-data_echo_12).*dt_echo_12/T_evaluated_echo_12#+
			#sum(abs2,_population_T1_01[:,1:3]-data_t1_01).*dt_t1_01/T_evaluated_t1_01+
			#sum(abs2,_population_T1_12[:,1:3]-data_t1_12).*dt_t1_12/T_evaluated_t1_12

	global _function_call
	_function_call += 1
	if(_function_call%5==0)
		println("Function call: ",_function_call," Loss = ",_loss)
	end
	return _loss
end

# initial guess
p_initial = [4.10817777;3.88303940; # transition frequencies
			 gamma1[1];gamma1[2]; # parameters for the decoupling
			 gamma2[1];gamma2[2]; # parameters for the dephasing
			 -1.017e-5;1.058e-4] # charge noise
#loss_gala(p_initial,[])
# upper and lower bound for optimizaiton
perturbation = [0.02; 0.02;  # freqs
		 		 0.9; 0.9;   # gam1
		 	 	 0.9; 0.9;	   # gma2
				-2.0; 1.0]   # charge noise

lower_bound = (1.0.-perturbation).*p_initial
upper_bound = (1.0.+perturbation).*p_initial
#@time Zygote.gradient(loss, p_initial)

########################################################
# Optimization
########################################################

optimization_function_az = OptimizationFunction(loss_gala, GalacticOptim.AutoZygote())
prob_az = GalacticOptim.OptimizationProblem(optimization_function_az, p_initial,
										    lb = lower_bound, ub = upper_bound)

_function_call = 0
_iteration_number = 0

@time sol = solve(prob_az , Opt(:LD_LBFGS,length(p_initial)),
				  maxiters=50,
				  ftol_rel=1e-3)

########################################################
# Visualization
########################################################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_echo_01,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 0,
				 TC_01,t_echo_01,N_states)
optim_echo_01 = GLOQ.get_population(rho_Echo01_u)
optim_echo_01 = optim_echo_01[:,1:3]

# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_echo_12,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 1,
				 TC_12,t_echo_12,N_states)
optim_echo_12 = GLOQ.get_population(rho_Echo12_u)
optim_echo_12 = optim_echo_12[:,1:3]

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_ramsey_01,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 0,
				 TC_01,t_ramsey_01,N_states)
optim_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)
optim_ramsey_01 = optim_ramsey_01[:,1:3]


# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_ramsey_12,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 1,
				 TC_12,t_ramsey_12,N_states)
optim_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)
optim_ramsey_12 = optim_ramsey_12[:,1:3]


# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_t1_01,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 0,
				 TC_01,t_t1_01,N_states)
optim_t1_01 = GLOQ.get_population(rho_T101_u)
optim_t1_01 = optim_t1_01[:,1:3]

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 [2.0*pi.*sol.u[1:2];omegas[3]],omr_t1_12,
				 [2.0*pi.*sol.u[7:8];char_noise[3]],
				 [sol.u[3:4];0.0],[sol.u[5:6];0.0],
				 1,
				 TC_12,t_t1_12,N_states)
optim_t1_12 = GLOQ.get_population(rho_T112_u)
optim_t1_12 = optim_t1_12[:,1:3]

# Plotting
fnt = Plots.font("Helvetica",16)
lfnt = Plots.font("Helvetica",12)
Plots.default(titlefont=fnt,
			  guidefont=fnt,
			  tickfont=fnt,
			  legendfont=lfnt,
			  linewidth=2);

# Echo 0-1
fig_e01 = plot(t_echo_01./1000.0,data_echo_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_e01,t_echo_01./1000.0,optim_echo_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_e01)

fig_e01_error = plot(t_ramsey_01./1000.0,optim_echo_01-data_echo_01,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_e01_error)

# Echo 1-2
fig_e12 = plot(t_echo_12./1000.0,data_echo_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_e12,t_echo_12./1000.0,optim_echo_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_e12)

fig_e12_error = plot(t_echo_12./1000.0,optim_echo_12-data_echo_12,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_e12_error)

# Ramsey 0-1
fig_r01 = plot(t_ramsey_01./1000.0,data_ramsey_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_r01,t_ramsey_01./1000.0,optim_ramsey_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_r01)

fig_r01_error = plot(t_ramsey_01./1000.0,optim_ramsey_01-data_ramsey_01,
		 size=(1200,600),legend=:outerright,
		 title = "Ramsey 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_r01_error)

# Ramsey 1-2
fig_r12 = plot(t_ramsey_12./1000.0,data_ramsey_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_r12,t_ramsey_12./1000.0,optim_ramsey_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_r12)

fig_r12_error = plot(t_ramsey_12./1000.0,optim_ramsey_12-data_ramsey_12,
				   size=(1200,600),legend=:outerright,
		   		   title="Ramsey 1-2",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_r12_error)

# T1 0-1
fig_t01 = plot(t_t1_01./1000.0,data_t1_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_t01,t_t1_01./1000.0,optim_t1_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_t01)

fig_t01_error = plot(t_t1_01./1000.0,optim_t1_01-data_t1_01,
		 size=(1200,600),legend=:outerright,
		 title = "T1 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_t01_error)

# T1 1-2
fig_t12 = plot(t_t1_12./1000.0,data_t1_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data"])
plot!(fig_t12,t_t1_12./1000.0,optim_t1_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_t12)

fig_t12_error = plot(t_t1_12./1000.0,optim_t1_12-data_t1_12,
		 size=(1200,600),legend=:outerright,
		 title = "T1 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$"])
display(fig_t12_error)

# output
println("Optimized results:",sol.u)
