using GLOQ, Zygote, LinearAlgebra
using GalacticOptim,NLopt,Optim
using Random
using Plots

Random.seed!(14);

# System parameters for a simple two level open quantum system
N_states = 2; # number of states
freqs = [4.1] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [1.0/(55.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
gamma2   = [1.0/(15.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)
omr = 2.0*pi*(4.1 - 5.0e-4) # drive frequency
TC = 2.5*17.0 # total control time in ns needed for a Rabi pulse 

# Initial state
initial_state = 0
state_u0 = [0.0;0.0]
state_v0 = [0.0;0.0]
state_u0[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time given in Microseconds
T_Ramsey = 40.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_times = 401
t_dark_times = collect(range(0.0, T_Ramsey, length=N_dark_times))

# Forward solve to generate synthetic data
rho_synthetic_ramsey_u,rho_synthetic_ramsey_v = GLOQ.RamseyForwardSolve(
				 state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
			     omegas,omr, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing parameters 
				 initial_state, # initial state
				 TC,t_dark_times,N_states) # control time, dark time, total number of states
population_synthetic = GLOQ.get_population(rho_synthetic_ramsey_u)

# Add noise to the synthetic data
noisy_data = copy(population_synthetic)
additive_noise = 0.025*randn(N_dark_times) 
noisy_data .+= additive_noise

# Physically, noisy data of population must be a probability.
# Shift and rescale the data so it is between [0,1]
for j = 1:N_states
	shift = minimum(noisy_data[:,j])
	if(shift<0.0)
		noisy_data[:,j] .-= shift
	end
end

for i = 1:N_dark_times
	noisy_data[i,:] ./= sum(noisy_data[i,:])
end

# Plot the synthetic data
fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
scatter!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,noisy_data)
display(fig)

# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#	 p[1] = transition frequency in GHz
#	 p[2] = gamma2
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p,dummy_parameter)
	_rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
				     (2*pi).*[p[1]],omr,
					 gamma1,[p[2]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_times,N_states)
	_population_ramsey = GLOQ.get_population(_rho_ramsey_u)

	_loss = sum(abs2,_population_ramsey-noisy_data)/N_dark_times
	return _loss
end

plot_callback = function(p,other_args)
	rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
					 (2*pi).*[p[1]],omr,
					 gamma1,[p[2]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_times,N_states)
	population_ramsey = GLOQ.get_population(rho_ramsey_u)
	fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,noisy_data,label=["Noisy-0" "Noisy-1"],
			 line = (:dash,0.0), marker = ([:hex :hex], 5, 0.5)  )
	plot!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_ramsey,label=["Opt-0" "Opt-1"]
		 )			#,line = (:dot, 4), marker = ([:hex :hex], 5, 0.1))
	display(fig)
    return false
end


p_true = [freqs;gamma2] # values to generate synthetic data
# initial guess for the optimization
p_initial = [freqs .- 0.5e-4;0.9.*gamma2]

# bounds for the optimization
lower_bound = (0.5).*p_true
upper_bound = (1.5).*p_true


# construct optimization object, use Zygote auto-differentiation to compute the gradient
loss_gradient = GalacticOptim.OptimizationFunction(loss, GalacticOptim.AutoZygote())
opt_prob = GalacticOptim.OptimizationProblem(loss_gradient, p_initial,
										 lb = lower_bound, ub = upper_bound)



println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = GalacticOptim.solve(opt_prob ,Fminbox(LBFGS()),
								cb = plot_callback,
								outer_iterations = 20,
								iterations = 10,
								show_trace=true,
								f_tol = 1e-10,
								outer_f_tol = 1e-10)
println("Optim Fminbox(LBFGS) Optimization done")
# present the solutions
println("\nOptimized results: ",sol.u,
        "\nLoss: ",sol.minimum,
		"\nError: ",sol.u-p_true)


println("NLopt LBFGS Optimization starts")
@time sol_nlopt_LBFGS = GalacticOptim.solve(opt_prob, Opt(:LD_LBFGS,length(p_initial)),
							  maxiters=200,
							  cb = plot_callback,
							  ftol_rel=1e-7)
println("NLopt LBFGS Optimization done")


# present the solution
println("\n\n\n---------------------------\nOptimized results: \n[Transition freq. (GHz), gamma2]\n",sol.u,
        "\nLoss: ",sol.minimum,
		"\nError: ",sol.u-p_true,"\n---------------------------\n")
