using Zygote
using LinearAlgebra
using NLopt
using Plots
using GLOQ

using Random
Random.seed!(14);

# System parameters for a simple two level open quantum system
N_states = 2; # number of states
freqs = [4.1] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [1.0/(45.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
gamma2   = [1.0/(24.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)
omr_ramsey = 2.0*pi*(4.1 - 5.0e-4) # drive frequency for the Ramsey experiment
omr_echo = 2.0*pi*4.1  # drive frequency for the Echo experiment
omr_t1 = 2.0*pi*4.1  # drive frequency for the T1 experiment
TC = 2.5*17.0 # total control time

# Initial state
initial_state = 0
state_u0 = [0.0;0.0]
state_v0 = [0.0;0.0]
state_u0[initial_state+1] = 1.0

# Forward solve to generate synthetic data for the Ramsey experiment
# Duration of the Ramsey experiment, largest dark time
T_Ramsey = 10.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_ramsey = 201
dt_ramsey = T_Ramsey/N_dark_ramsey
t_dark_ramsey = collect(range(0.0, T_Ramsey, length=N_dark_ramsey))
# Forward solve
rho_synthetic_ramsey_u,rho_synthetic_ramsey_v = GLOQ.RamseyForwardSolve(
				 state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
			     omegas,omr_ramsey, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing 
				 initial_state, # initial state
				 TC,t_dark_ramsey,N_states) # control time, dark time, total number of states
population_ramsey_synthetic = GLOQ.get_population(rho_synthetic_ramsey_u)

# Echo experiment
# Duration of the Echo experiment, largest dark time
T_Echo = 30.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_echo = 601
# Step size of the Echo experiment
dt_echo = T_Echo/N_dark_echo
t_dark_echo = collect(range(0.0, T_Echo, length=N_dark_echo))
rho_synthetic_echo_u,rho_synthetic_echo_v = GLOQ.EchoForwardSolve(
				 state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
			     omegas,omr_echo, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing 
				 initial_state, # initial state
				 TC,t_dark_echo,N_states) # control time, dark time, total number of states
population_echo_synthetic = GLOQ.get_population(rho_synthetic_echo_u)

# T1
# Duration of the T1-decay experiment, largest dark time
T_t1 = 40.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_t1 = 401
# Step size of the T1 experiment
dt_t1 = T_t1/N_dark_t1
t_dark_t1 = collect(range(0.0, T_t1, length=N_dark_t1))
rho_synthetic_t1_u,rho_synthetic_t1_v = GLOQ.T1ForwardSolve(
				 state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
			     omegas,omr_t1, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing 
				 initial_state, # initial state
				 TC,t_dark_t1,N_states) # control time, dark time, total number of states
population_t1_synthetic = GLOQ.get_population(rho_synthetic_t1_u)


# Plot the synthetic data
#fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
#display(fig)

# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#	 p[1] = transition frequency in GHz
#    p[2] = gamma1
#    p[3] = gamma2
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p)
	# Ramsey
	_rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
				     (2*pi).*[p[1]],omr_ramsey,
					 [p[2]],[p[3]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_ramsey,N_states)
	_population_ramsey = GLOQ.get_population(_rho_ramsey_u)
	# Echo
    _rho_echo_u,_rho_echo_v = GLOQ.EchoForwardSolve(state_u0,state_v0,
				     (2*pi).*[p[1]],omr_echo,
					 [p[2]],[p[3]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_echo,N_states)
	_population_echo = GLOQ.get_population(_rho_echo_u)
	# T1
    _rho_t1_u,_rho_t1_v = GLOQ.T1ForwardSolve(state_u0,state_v0,
				     (2*pi).*[p[1]],omr_t1,
					 [p[2]],[p[3]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_t1,N_states)
	_population_t1 = GLOQ.get_population(_rho_t1_u)

	_loss = sum(abs2,_population_ramsey-population_ramsey_synthetic)*dt_ramsey+
			sum(abs2,_population_echo-population_echo_synthetic)*dt_echo+
			sum(abs2,_population_t1-population_t1_synthetic)*dt_t1
		
	return _loss
end

function loss_nlopt(p,grad)
	grad .= Zygote.gradient(loss,p)[1]
	return loss(p)
end


p_true = [freqs;gamma1;gamma2] # values to generate synthetic data
# initial guess for the optimization
p_initial = [freqs.-0.5e-4;0.75.*gamma1;1.25.*gamma2]
# bounds for the optimization
lower_bound = (0.5).*p_true
upper_bound = (1.5).*p_true


println("NLopt LBFGS Optimization starts")
nlopt_obj = NLopt.Opt(:LD_LBFGS,length(p_initial))
nlopt_obj.min_objective = loss_nlopt
nlopt_obj.lower_bounds = lower_bound
nlopt_obj.upper_bounds = upper_bound
nlopt_obj.maxeval = 200
nlopt_obj.ftol_rel = 1e-4
@time loss_value_nlopt,sol_nlopt,flag_nlopt = NLopt.optimize(nlopt_obj,p_initial)
println("NLopt LBFGS Optimization done")


# present the solution
println("\n\n\n---------------------------\nOptimized results: \n[Transition freq. (GHz), gamma2]\n",sol_nlopt,
        "\nLoss: ",loss_value_nlopt,
		"\nError: ",sol_nlopt-p_true,"\n---------------------------\n")
		