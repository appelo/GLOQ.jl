using Zygote
using Random
using DelimitedFiles
using LinearAlgebra
using GalacticOptim,Optim,NLopt
using DifferentialEquations

using Plots
using LaTeXStrings
#using GLOQ
#using ReverseDiff
include("../src/GLOQ.jl")
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

########################################################################
# Step 1: set up the problem
########################################################################
N_states = 4;
# True values to generate the synthtetic data
freqs_true = [4.0108; 3.8830; 3.6287];
omegas_true = 2.0*pi.*freqs_true;
gamma1_true   = [2.222222e-05; 4.761904e-05; 4.545455e-5];
gamma2_true   = [4.166667e-05; 6.8965517e-5; 2.3255814e-4];
char_noise_true = [2e-6,6e-5,1.3e-4];



########################################################################
# Step 1a: known parameters to set up the experiment
########################################################################
# Konwn driving frequencies
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

# Convert the dark time from micro-seconds to nano-seconds
T_T1_01 = 125.0*GLOQ.GLOQ_MICRO_SEC
T_T1_12 = 125.0*GLOQ.GLOQ_MICRO_SEC
T_T1_23 = 125.0*GLOQ.GLOQ_MICRO_SEC
T_Ramsey_01 = 20.0*GLOQ.GLOQ_MICRO_SEC
T_Ramsey_12 = 20.0*GLOQ.GLOQ_MICRO_SEC
T_Ramsey_23 = 10.0*GLOQ.GLOQ_MICRO_SEC
T_Echo_01 = 40.0*GLOQ.GLOQ_MICRO_SEC
T_Echo_12 = 20.0*GLOQ.GLOQ_MICRO_SEC
T_Echo_23 = 10.0*GLOQ.GLOQ_MICRO_SEC

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


# Knwon initial conditions for the quantum system
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

########################################################################
# Step 1b: read the synthetic data
########################################################################
data_t1_01 = readdlm("synthetic_data/t1_01.txt",',')
data_t1_12 = readdlm("synthetic_data/t1_12.txt",',')
data_t1_23 = readdlm("synthetic_data/t1_23.txt",',')
data_ramsey_01 = readdlm("synthetic_data/ramsey_01.txt",',')
data_ramsey_12 = readdlm("synthetic_data/ramsey_12.txt",',')
data_ramsey_23 = readdlm("synthetic_data/ramsey_23.txt",',')
data_echo_01 = readdlm("synthetic_data/echo_01.txt",',')
data_echo_12 = readdlm("synthetic_data/echo_12.txt",',')
data_echo_23 = readdlm("synthetic_data/echo_23.txt",',')

########################################################################
# add multiplicative and additive noise to the synthetic data
########################################################################
multiplicative_noise = 1.0.+0.025*randn(n1use)
additive_noise = 0.02*randn(n1use)
data_t1_01 .*= multiplicative_noise
data_t1_01 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n2use)
additive_noise = 0.025*randn(n2use)
data_t1_12 .*= multiplicative_noise
data_t1_12 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n3use)
additive_noise = 0.025*randn(n3use)
data_t1_23 .*= multiplicative_noise
data_t1_23 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n4use)
additive_noise = 0.025*randn(n4use)
data_ramsey_01 .*= multiplicative_noise
data_ramsey_01 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n5use)
additive_noise = 0.025*randn(n5use)
data_ramsey_12 .*= multiplicative_noise
data_ramsey_12 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n6use)
additive_noise = 0.025*randn(n6use)
data_ramsey_23 .*= multiplicative_noise
data_ramsey_23 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n7use)
additive_noise = 0.025*randn(n7use)
data_echo_01 .*= multiplicative_noise
data_echo_01 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n8use)
additive_noise = 0.025*randn(n8use)
data_echo_12 .*= multiplicative_noise
data_echo_12 .+= additive_noise

multiplicative_noise = 1.0.+0.025*randn(n9use)
additive_noise = 0.025*randn(n9use)
data_echo_23 .*= multiplicative_noise
data_echo_23 .+= additive_noise
# population (probability) must be in the range of [0,1],
# and the summation of the population for different states must be 1
function correct_noisy_data!(data)
	_n,_m = size(data)
	for i = 1:_n
		for j = 1:_m
			if(data[i,j]<0.0)
				data[i,j] = 0.0
			end
		end
		data[i,:] ./= sum(data[i,:])
	end
end

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
# Step 2: define and solve the optimization problem
########################################################
# Step 2a: define the loss (objective) function
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


########################################################################
# Step 2b: initial guess and bounds for target parameters
########################################################################
# transition frequency
freqs = [4.0108; 3.8830; 3.6287];
omegas = 2.0*pi.*freqs # change things from GHz to radians
# magnitude of the charge noise
char_noise = [1.5e-6, 5.0e-5, 1.0e-4]
# Lindblad parameters
gamma1   = [1.0e-05; 4.0e-05; 4.0e-5]
gamma2   = [2.5e-05; 5.0e-05; 2.5e-4]
# initial guess
p_initial = [freqs; # transition frequencies
			 char_noise;
			 gamma1; # parameters for the decoupling
			 gamma2 # parameters for the dephasing
 		 	 ]# charge noise

# upper and lower bound for optimizaiton
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

#####################################################################################
# Step 2c: optimization interface
#####################################################################################
# interface for the GalacticOptim with gradient computed by auto-differentiation
# (with the Zygote package)
optimization_function_az = OptimizationFunction(loss_gala, GalacticOptim.AutoZygote())
prob_az = GalacticOptim.OptimizationProblem(optimization_function_az, p_initial;
										    lb = lower_bound, ub = upper_bound)
_function_call = 0
# solve the optimization problem
@time sol = solve(prob_az , Opt(:LD_LBFGS,length(p_initial)),
				  maxiters=200,
				  ftol_rel=1e-3
				 )
p_optim = sol.u

########################################################
# Step 3: visualization
########################################################
# Plot results
include("ShowResultsOfCharacterizationSynthetic.jl")
# output
println("Transition frequencies: ", "\n Optimzied: ",p_optim[1:3], "\n True: ",freqs_true,"\n Error: ",p_optim[1:3]-freqs_true,"\n")
println("Charge noise:\n ", " Optimzied: ",p_optim[4:6], "\n True: ",char_noise_true,"\n Error: ",p_optim[4:6]-char_noise_true,"\n")
println("Gamma 1:\n ","\n Optimzed: ",p_optim[7:9],"\n True: ",gamma1_true,"\n Error: ",p_optim[7:9]-gamma1_true,"\n")
println("Gamma 2:\n ","\n Optimzed: ",p_optim[10:12],"\n True: ",gamma2_true,"\n Error: ",p_optim[10:12]-gamma1_true,"\n")
