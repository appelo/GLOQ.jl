using GLOQ, Zygote, LinearAlgebra
using Optim,NLopt
using Plots

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

# Plot the synthetic data
fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
display(fig)

# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#	 p[1] = transition frequency in GHz
#	 p[2] = gamma2
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss_optim(p)
	_rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
				     (2*pi).*[p[1]],omr,
					 gamma1,[p[2]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_times,N_states)
	_population_ramsey = GLOQ.get_population(_rho_ramsey_u)

	_loss = sum(abs2,_population_ramsey-population_synthetic)/N_dark_times
	return _loss
end

function gradient_optim!(G,p)
	G .= Zygote.gradient(loss_optim,p)[1]
end

function loss_nlopt(p,grad)
	grad .= Zygote.gradient(loss_optim,p)[1]
	return loss_optim(p)
end

plot_callback = function(_trace)
	p = _trace.metadata["x"]
	rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
					 (2*pi).*[p[1]],omr,
					 gamma1,[p[2]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_times,N_states)
	population_ramsey = GLOQ.get_population(rho_ramsey_u)
	fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic,label=["Syn-0" "Syn-1"],
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
println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = Optim.optimize(loss_optim,gradient_optim!,
						   lower_bound,upper_bound,p_initial,
						   Fminbox(LBFGS()),
						   Optim.Options(
							callback = plot_callback,extended_trace=true,
							outer_iterations = 10,
							iterations = 10,
							show_trace=true,
							f_tol = 1e-4,
							outer_f_tol = 1e-4))
println("Optim Fminbox(LBFGS) Optimization done")
# present the solutions
println("\nOptimized results: ",sol.minimizer,
        "\nLoss: ",sol.minimum,
		"\nError: ",sol.minimizer-p_true)


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
