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
# Plotting
fnt = Plots.font("Helvetica",16)
lfnt = Plots.font("Helvetica",12)
Plots.default(titlefont=fnt,
			  guidefont=fnt,
			  tickfont=fnt,
			  legendfont=lfnt,
			  linewidth=2);
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

N_states = 4;
########################################################################
freqs_true = [4.0108; 3.8830; 3.6287];
omegas_true = 2.0*pi.*freqs_true;
gamma1_true   = [2.222222e-05; 4.761904e-05; 4.545455e-5];
gamma2_true   = [4.166667e-05; 6.8965517e-5; 2.3255814e-4];
char_noise_true = [2e-6,6e-5,1.3e-4];

########################################################################
#=
freqs = [4.0108; 3.8830; 3.6287]
omegas = 2.0*pi.*freqs
char_noise = [1.5e-6, 5.0e-5, 1.29e-4]
gamma1   = [2.25e-05; 4.5e-05; 4.25e-5]
gamma2   = [4.0e-05; 7.0e-05; 2.0e-4]
=#

freqs = [4.0108; 3.8830; 3.6287];
omegas = 2.0*pi.*freqs
char_noise = [1.5e-6, 5.0e-5, 1.0e-4]
gamma1   = [1.0e-05; 4.0e-05; 4.0e-5]
gamma2   = [2.5e-05; 5.0e-05; 2.5e-4]

# driving frequency
omr_t1_01 = 2.0*pi*4.0108
omr_t1_12 = 2.0*pi*3.8830
omr_t1_23 = 2.0*pi*3.6287

omr_ramsey_01 = 2.0*pi*(4.0108-5e-4)
omr_ramsey_12 = 2.0*pi*(3.8830-2.0e-3)
omr_ramsey_23 = 2.0*pi*(3.6287-2.0e-3)

omr_echo_01 = 2.0*pi*4.0108
omr_echo_12 = 2.0*pi*3.8830
omr_echo_23 = 2.0*pi*3.6287




# width of the control signal
width = 21.25
TC_01 = 2.5*width

width = 17.00
TC_12 = 2.5*width

width = 17.00
TC_23 = 2.5*width

# in nano second
micro_sec = 1000.0 # probably include it in our source code
T_T1_01 = 125.0*micro_sec
T_T1_12 = 125.0*micro_sec
T_T1_23 = 125.0*micro_sec
T_Ramsey_01 = 20.0*micro_sec
T_Ramsey_12 = 20.0*micro_sec
T_Ramsey_23 = 10.0*micro_sec
T_Echo_01 = 40.0*micro_sec
T_Echo_12 = 20.0*micro_sec
T_Echo_23 = 10.0*micro_sec

# data size used for characterization
n1use = 201
n2use = 201
n3use = 201
n4yse = 201
n4use = 201
n5use = 201
n6use = 201
n7use = 201
n8use = 201
n9use = 201

# time points evaluated
dt_t1_01 = T_T1_01/200
dt_t1_12 = T_T1_12/200
dt_t1_23 = T_T1_23/200
dt_ramsey_01 = T_Ramsey_01/200
dt_ramsey_12 = T_Ramsey_12/200
dt_ramsey_23= T_Ramsey_23/200
dt_echo_01 = T_Echo_01/200
dt_echo_12 = T_Echo_12/200
dt_echo_23 = T_Echo_23/200

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

t_t1_23 = zeros(Float64,n3use)
for i = 1:n3use
	t_t1_23[i] = (i-1)*dt_t1_23
end
T_evaluated_t1_23 = t_t1_23[end]

t_ramsey_01 = zeros(Float64,n4use)
for i = 1:n4use
	t_ramsey_01[i] = (i-1)*dt_ramsey_01
end
T_evaluated_ramsey_01 = t_ramsey_01[end]

t_ramsey_12 = zeros(Float64,n5use)
for i = 1:n5use
	t_ramsey_12[i] = (i-1)*dt_ramsey_12
end
T_evaluated_ramsey_12 = t_ramsey_12[end]

t_ramsey_23 = zeros(Float64,n6use)
for i = 1:n6use
	t_ramsey_23[i] = (i-1)*dt_ramsey_23
end
T_evaluated_ramsey_23 = t_ramsey_23[end]

t_echo_01 = zeros(Float64,n7use)
for i = 1:n7use
	t_echo_01[i] = (i-1)*dt_echo_01
end
T_evaluated_echo_01 = t_echo_01[end]

t_echo_12 = zeros(Float64,n8use)
for i = 1:n8use
	t_echo_12[i] = (i-1)*dt_echo_12
end
T_evaluated_echo_12 = t_echo_12[end]

t_echo_23 = zeros(Float64,n9use)
for i = 1:n9use
	t_echo_23[i] = (i-1)*dt_echo_23
end
T_evaluated_echo_23 = t_echo_23[end]


# initial conditions
u0_t1_01 = [1.0;0.0;0.0;0.0]
u0_t1_12 = [0.0;1.0;0.0;0.0]
u0_t1_23 = [0.0;0.0;1.0;0.0]
u0_ramsey_01 = [1.0;0.0;0.0;0.0]
u0_ramsey_12 = [0.0;1.0;0.0;0.0]
u0_ramsey_23 = [0.0;0.0;1.0;0.0]
u0_echo_01 = [1.0;0.0;0.0;0.0]
u0_echo_12 = [0.0;1.0;0.0;0.0]
u0_echo_23 = [0.0;0.0;1.0;0.0]

v0_t1_01 = [0.0;0.0;0.0;0.0]
v0_t1_12 = [0.0;0.0;0.0;0.0]
v0_t1_23 = [0.0;0.0;0.0;0.0]
v0_ramsey_01 = [0.0;0.0;0.0;0.0]
v0_ramsey_12 = [0.0;0.0;0.0;0.0]
v0_ramsey_23 = [0.0;0.0;0.0;0.0]
v0_echo_01 = [0.0;0.0;0.0;0.0]
v0_echo_12 = [0.0;0.0;0.0;0.0]
v0_echo_23 = [0.0;0.0;0.0;0.0]

data_t1_01 = readdlm("synthetic_data/t1_01.txt",',')
data_t1_12 = readdlm("synthetic_data/t1_12.txt",',')
data_t1_23 = readdlm("synthetic_data/t1_23.txt",',')
data_ramsey_01 = readdlm("synthetic_data/ramsey_01.txt",',')
data_ramsey_12 = readdlm("synthetic_data/ramsey_12.txt",',')
data_ramsey_23 = readdlm("synthetic_data/ramsey_23.txt",',')
data_echo_01 = readdlm("synthetic_data/echo_01.txt",',')
data_echo_12 = readdlm("synthetic_data/echo_12.txt",',')
data_echo_23 = readdlm("synthetic_data/echo_23.txt",',')

# add noise to the synthetic data
multiplicative_noise = 1.0.+0.025*randn(201)
data_t1_01 .*= multiplicative_noise
data_t1_12 .*= multiplicative_noise
data_t1_23 .*= multiplicative_noise
data_ramsey_01 .*= multiplicative_noise
data_ramsey_12 .*= multiplicative_noise
data_ramsey_23 .*= multiplicative_noise
data_echo_01 .*= multiplicative_noise
data_echo_12 .*= multiplicative_noise
data_echo_23 .*= multiplicative_noise

# population must be in the range of [0,1]
function correct_noisy_data!(data)
	_n,_m = size(data)
	for i = 1:_n
		for j = 1:_m
			if(data[i,j]>1.0)
				data[i,j] = 1.0
			elseif(data[i,j]<0.0)
				data[i,j] = 0.0
			end
		end
	end
end
# normalize the data

correct_noisy_data!(data_t1_01)
correct_noisy_data!(data_t1_12)
correct_noisy_data!(data_t1_23)
correct_noisy_data!(data_echo_01)
correct_noisy_data!(data_echo_12)
correct_noisy_data!(data_echo_23)
correct_noisy_data!(data_ramsey_01)
correct_noisy_data!(data_ramsey_12)
correct_noisy_data!(data_ramsey_23)

global _function_call = 0
global _iteration_number = 0
########################################################
# Define loss function
########################################################
function loss_gala(p,p_keywords)
	# Ramsey experiments
	_rho_Ramsey_01_u,_rho_Ramsey_01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
					 2.0*pi.*p[1:3],omr_ramsey_01,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 0,
					 TC_01,t_ramsey_01,N_states)
	_population_Ramsey_01 = GLOQ.get_population(_rho_Ramsey_01_u)

	_rho_Ramsey_12_u,_rho_Ramsey_12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
					 2.0*pi.*p[1:3],omr_ramsey_12,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 1,
					 TC_12,t_ramsey_12,N_states)
	_population_Ramsey_12 = GLOQ.get_population(_rho_Ramsey_12_u)

	_rho_Ramsey_23_u,_rho_Ramsey_23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
					 2.0*pi.*p[1:3],omr_ramsey_23,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 2,
					 TC_23,t_ramsey_23,N_states)
	_population_Ramsey_23 = GLOQ.get_population(_rho_Ramsey_23_u)


	# Echo experiments
	_rho_Echo_01_u,_rho_Echo_01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
					 2.0*pi.*p[1:3],omr_echo_01,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 0,
					 TC_01,t_echo_01,N_states)
	_population_Echo_01 = GLOQ.get_population(_rho_Echo_01_u)

	_rho_Echo_12_u,_rho_Echo_12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
					 2.0*pi.*p[1:3],omr_echo_12,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 1,
					 TC_12,t_echo_12,N_states)
	_population_Echo_12 = GLOQ.get_population(_rho_Echo_12_u)

	_rho_Echo_23_u,_rho_Echo_23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
					 2.0*pi.*p[1:3],omr_echo_23,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 2,
					 TC_23,t_echo_23,N_states)
	_population_Echo_23 = GLOQ.get_population(_rho_Echo_23_u)

	# T1-decay experiments
	_rho_T1_01_u,_rho_T1_01_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
					 2.0*pi.*p[1:3],omr_t1_01,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 0,
					 TC_01,t_t1_01,N_states)
	_population_T1_01 = GLOQ.get_population(_rho_T1_01_u)

	_rho_T1_12_u,_rho_T1_12_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
					 2.0*pi.*p[1:3],omr_t1_12,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 1,
					 TC_12,t_t1_12,N_states)
	_population_T1_12 = GLOQ.get_population(_rho_T1_12_u)

	_rho_T1_23_u,_rho_T1_23_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
					 2.0*pi.*p[1:3],omr_t1_23,
					 2.0*pi.*p[4:6],
					 p[7:9],p[10:12],
					 2,
					 TC_23,t_t1_23,N_states)
	_population_T1_23 = GLOQ.get_population(_rho_T1_23_u)

	#################################################################################
#	_loss = sum(abs2,_population_Echo_23-data_echo_23).*dt_echo_23/T_evaluated_echo_23

	_loss = sum(abs2,_population_Ramsey_01-data_ramsey_01)*dt_ramsey_01/T_evaluated_ramsey_01+
			sum(abs2,_population_Ramsey_12-data_ramsey_12)*dt_ramsey_12/T_evaluated_ramsey_12+
			sum(abs2,_population_Ramsey_23-data_ramsey_23)*dt_ramsey_23/T_evaluated_ramsey_23+
			sum(abs2,_population_Echo_01-data_echo_01).*dt_echo_01/T_evaluated_echo_01+
			sum(abs2,_population_Echo_12-data_echo_12).*dt_echo_12/T_evaluated_echo_12+
			sum(abs2,_population_Echo_23-data_echo_23).*dt_echo_23/T_evaluated_echo_23+
			sum(abs2,_population_T1_01-data_t1_01).*dt_t1_01/T_evaluated_t1_01+
			sum(abs2,_population_T1_12-data_t1_12).*dt_t1_12/T_evaluated_t1_12+
			sum(abs2,_population_T1_23-data_t1_23).*dt_t1_23/T_evaluated_t1_23

	global _function_call
	_function_call += 1
	if(_function_call%1==0)
		println("Function call: ",_function_call," Loss = ",_loss)
	end
	return _loss
end

########################################################
# Optimization
########################################################
# initial guess
p_initial = [freqs; # transition frequencies
			 char_noise;
			 gamma1; # parameters for the decoupling
			 gamma2 # parameters for the dephasing
 		 	 ]# charge noise
####

# upper and lower bound for optimizaiton

#=
p_true = [freqs_true; # transition frequencies
			char_noise_true;
			 gamma1_true; # parameters for the decoupling
			 gamma2_true # parameters for the dephasing
 		 	 ]
p_initial = p_true#0.999.*p_true
=#
initial_loss = loss_gala(p_initial,[])

lower_bound = [ 4.00; 3.85; 3.60;
				1e-6; 5e-5; 1e-4;
				1e-5; 3e-5; 3e-5;
				2e-5; 2.5e-5; 2e-4;
			  ]
upper_bound = [ 4.10; 3.95; 3.65;
				3e-6; 7e-5; 2e-4;
				3e-5; 5e-5; 5e-5;
				5e-5; 1e-4; 3e-4;
			  ]


optimization_function_az = OptimizationFunction(loss_gala, GalacticOptim.AutoZygote())
prob_az = GalacticOptim.OptimizationProblem(optimization_function_az, p_initial;#p_12,#p_initial,#p_12,
										    lb = lower_bound, ub = upper_bound)

_function_call = 0
_iteration_number = 0

@time sol = solve(prob_az , Opt(:LD_LBFGS,length(p_initial)),
				  maxiters=200,
				  ftol_rel=1e-3)


p_optim = sol.u

########################################################
# Visualization
########################################################

########################
# Initial parameters
########################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 2.0*pi.*p_initial[1:3],omr_echo_01,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 0,
				 TC_01,t_echo_01,N_states)
initial_echo_01 = GLOQ.get_population(rho_Echo01_u)


# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 2.0*pi.*p_initial[1:3],omr_echo_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_echo_12,N_states)
initial_echo_12 = GLOQ.get_population(rho_Echo12_u)

# Echo 2-3
rho_Echo23_u,rho_Echo23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
				 2.0*pi.*p_initial[1:3],omr_echo_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_echo_23,N_states)
initial_echo_23 = GLOQ.get_population(rho_Echo23_u)

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 2.0*pi.*p_optim[1:3],omr_ramsey_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_ramsey_01,N_states)
initial_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)

# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 2.0*pi.*p_initial[1:3],omr_ramsey_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_ramsey_12,N_states)
initial_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)

# Ramsey 2-3
rho_Ramsey23_u,rho_Ramsey23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
				 2.0*pi.*p_initial[1:3],omr_ramsey_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_ramsey_23,N_states)
initial_ramsey_23 = GLOQ.get_population(rho_Ramsey23_u)


# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 2.0*pi.*p_initial[1:3],omr_t1_01,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 0,
				 TC_01,t_t1_01,N_states)
initial_t1_01 = GLOQ.get_population(rho_T101_u)

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 2.0*pi.*p_initial[1:3],omr_t1_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_t1_12,N_states)
initial_t1_12 = GLOQ.get_population(rho_T112_u)

# T1 2-3
rho_T123_u,rho_T123_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
				 2.0*pi.*p_initial[1:3],omr_t1_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_t1_23,N_states)
initial_t1_23 = GLOQ.get_population(rho_T123_u)

# Echo 0-1
initial_e01_error = plot(t_ramsey_01./1000.0,initial_echo_01-data_echo_01,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e01_error)

# Echo 1-2
initial_e12_error = plot(t_echo_12./1000.0,initial_echo_12-data_echo_12,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e12_error)

# Echo 2-3
initial_e23 = plot(t_echo_23./1000.0,data_echo_23,
		 size=(1200,600),legend=:outerright,line=(:dash),
		 title = "Echo 2-3",
		 label=[L"$\rho_{00}$ data" L"$\rho_{11}$ data" L"$\rho_{22}$ data" L"$\rho_{33}$ data"])
plot!(initial_e23,t_echo_23./1000.0,initial_echo_23,
	  size=(1200,600),legend=:outerright,
	  title = "Echo 2-3",
	  label=[L"$\rho_{00}$ sim" L"$\rho_{11}$ sim" L"$\rho_{22}$ sim" L"$\rho_{33}$ sim"])
display(initial_e23)

initial_e23_error = plot(t_echo_23./1000.0,initial_echo_23-data_echo_23,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e23_error)

# Ramsey 0-1
initial_r01_error = plot(t_ramsey_01./1000.0,initial_ramsey_01-data_ramsey_01,
		 size=(1200,600),legend=:outerright,
		 title = "Ramsey 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r01_error)

# Ramsey 1-2
initial_r12_error = plot(t_ramsey_12./1000.0,initial_ramsey_12-data_ramsey_12,
				   size=(1200,600),legend=:outerright,
		   		   title="Ramsey 1-2",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r12_error)

# Ramsey 2-3
initial_r23_error = plot(t_ramsey_23./1000.0,initial_ramsey_23-data_ramsey_23,
				   size=(1200,600),legend=:outerright,
		   		   title="Ramsey 2-3",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r23_error)

# T1 0-1
initial_t01_error = plot(t_t1_01./1000.0,initial_t1_01-data_t1_01,
		 size=(1200,600),legend=:outerright,
		 title = "T1 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t01_error)

# T1 1-2
initial_t12_error = plot(t_t1_12./1000.0,initial_t1_12-data_t1_12,
		 size=(1200,600),legend=:outerright,
		 title = "T1 1-2",
		 		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t12_error)

# T1 2-3
initial_t23_error = plot(t_t1_12./1000.0,initial_t1_23-data_t1_23,
		 size=(1200,600),legend=:outerright,
		 title = "T1 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t23_error)

########################
# Optimized parameters
########################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 2.0*pi.*p_optim[1:3],omr_echo_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_echo_01,N_states)
optim_echo_01 = GLOQ.get_population(rho_Echo01_u)


# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 2.0*pi.*p_optim[1:3],omr_echo_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_echo_12,N_states)
optim_echo_12 = GLOQ.get_population(rho_Echo12_u)

# Echo 2-3
rho_Echo23_u,rho_Echo23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
				 2.0*pi.*p_optim[1:3],omr_echo_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_echo_23,N_states)
optim_echo_23 = GLOQ.get_population(rho_Echo23_u)

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 2.0*pi.*p_optim[1:3],omr_ramsey_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_ramsey_01,N_states)
optim_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)

# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 2.0*pi.*p_optim[1:3],omr_ramsey_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_ramsey_12,N_states)
optim_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)

# Ramsey 2-3
rho_Ramsey23_u,rho_Ramsey23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
				 2.0*pi.*p_optim[1:3],omr_ramsey_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_ramsey_23,N_states)
optim_ramsey_23 = GLOQ.get_population(rho_Ramsey23_u)


# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 2.0*pi.*p_optim[1:3],omr_t1_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_t1_01,N_states)
optim_t1_01 = GLOQ.get_population(rho_T101_u)

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 2.0*pi.*p_optim[1:3],omr_t1_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_t1_12,N_states)
optim_t1_12 = GLOQ.get_population(rho_T112_u)

# T1 2-3
rho_T123_u,rho_T123_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
				 2.0*pi.*p_optim[1:3],omr_t1_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_t1_23,N_states)
optim_t1_23 = GLOQ.get_population(rho_T123_u)



# Echo 0-1
fig_e01 = plot(t_echo_01./1000.0,data_echo_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_e01,t_echo_01./1000.0,optim_echo_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"])
display(fig_e01)

fig_e01_error = plot(t_ramsey_01./1000.0,optim_echo_01-data_echo_01,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_e01_error)

# Echo 1-2
fig_e12 = plot(t_echo_12./1000.0,data_echo_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_e12,t_echo_12./1000.0,optim_echo_12,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_e12)

fig_e12_error = plot(t_echo_12./1000.0,optim_echo_12-data_echo_12,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_e12_error)

# Echo 2-3
fig_e23 = plot(t_echo_23./1000.0,data_echo_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_e23,t_echo_23./1000.0,optim_echo_23,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_e23)

fig_e23_error = plot(t_echo_23./1000.0,optim_echo_23-data_echo_23,
		 size=(1200,600),legend=:outerright,
		 title = "Echo 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_e23_error)

# Ramsey 0-1
fig_r01 = plot(t_ramsey_01./1000.0,data_ramsey_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_r01,t_ramsey_01./1000.0,optim_ramsey_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r01)

fig_r01_error = plot(t_ramsey_01./1000.0,optim_ramsey_01-data_ramsey_01,
		 size=(1200,600),legend=:outerright,
		 title = "Ramsey 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_r01_error)

# Ramsey 1-2
fig_r12 = plot(t_ramsey_12./1000.0,data_ramsey_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_r12,t_ramsey_12./1000.0,optim_ramsey_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r12)

fig_r12_error = plot(t_ramsey_12./1000.0,optim_ramsey_12-data_ramsey_12,
				   size=(1200,600),legend=:outerright,
		   		   title="Ramsey 1-2",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_r12_error)

# Ramsey 2-3
fig_r23 = plot(t_ramsey_23./1000.0,data_ramsey_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_r23,t_ramsey_23./1000.0,optim_ramsey_23,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r23)

fig_r23_error = plot(t_ramsey_23./1000.0,optim_ramsey_23-data_ramsey_23,
				   size=(1200,600),legend=:outerright,
		   		   title="Ramsey 2-3",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_r23_error)

# T1 0-1
fig_t01 = plot(t_t1_01./1000.0,data_t1_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_t01,t_t1_01./1000.0,optim_t1_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t01)

fig_t01_error = plot(t_t1_01./1000.0,optim_t1_01-data_t1_01,
		 size=(1200,600),legend=:outerright,
		 title = "T1 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_t01_error)

# T1 1-2
fig_t12 = plot(t_t1_12./1000.0,data_t1_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_t12,t_t1_12./1000.0,optim_t1_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t12)

fig_t12_error = plot(t_t1_12./1000.0,optim_t1_12-data_t1_12,
		 size=(1200,600),legend=:outerright,
		 title = "T1 1-2",
		 		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_t12_error)

# T1 2-3
fig_t23 = plot(t_t1_23./1000.0,data_t1_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_t23,t_t1_23./1000.0,optim_t1_23,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t23)

fig_t23_error = plot(t_t1_23./1000.0,optim_t1_23-data_t1_23,
		 size=(1200,600),legend=:outerright,
		 title = "T1 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_t23_error)

# output
println("Transition frequencies: ", " optimzied: ",p_optim[1:3], " true: ",freqs_true)
println("Charge noise: ", " optimzied: ",p_optim[4:6], " true: ",char_noise_true)
println("Gamma 1: "," optimzed: ",p_optim[7:9]," true: ",gamma1_true)
println("Gamma 2: "," optimzed: ",p_optim[10:12]," true: ",gamma2_true)
