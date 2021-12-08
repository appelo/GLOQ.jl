
## Example 1: characterization of a single qubit Ramsey experiment without noise
In this example, we solve a characterization problem for a single qubit Ramsey experiment. We seek the dephasing time $T_2=1/\gamma_2$ and the transition frequency $\omega_{01}$.
This example is corresponding to `GLOQ.jl/examples/SingleQubit.jl`. 
The code can be invoked by e.g. `cd("examples");include("SingleQubit.jl")`. 
### Step 1: generate the synthetic data
```julia
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
rho_u0 = [0.0;0.0]
rho_v0 = [0.0;0.0]
rho_u0[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time given in Microseconds
T_Ramsey = 40.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_times = 401
t_dark_times = collect(range(0.0, T_Ramsey, length=N_dark_times))

# Forward solve to generate synthetic data
rho_synthetic_ramsey_u,rho_synthetic_ramsey_v = GLOQ.RamseyForwardSolve(
				 rho_u0,rho_v0, # initial values, u for the real part, v for the imaginary part
				 omegas,omr, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing parameters 
				 initial_state, # initial state
				 TC,t_dark_times,N_states) # control time, dark time, total number of states
population_synthetic = GLOQ.get_population(rho_synthetic_ramsey_u)
```
### Step 2: define the objective function, initial guess and optimization bounds
#### Step 2a: define the cost objective function. Here, we use the normalized l2-mismatch. 

Suppose the number of dark times in the Ramsey experiment to be $N_{\textrm{Dark}}$. The objective function is defined as 
```math
\big|\big| \textrm{Forward Solve Results}-\textrm{Synthetic Data}\big|\big|_2^2/N_{\textrm{Dark}},
```
with $||\cdot||_2$ being the standard $l_2$ norm.
```julia
# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#	 p[1] = transition frequency in GHz
#	 p[2] = gamma2
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p,dummy_parameter)
	_rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
					 (2*pi).*[p[1]],omr,
					 gamma1,[p[2]],#gamma1,gamma2,
					 initial_state, # initial state
					 TC,t_dark_times,N_states)
	_population_ramsey = GLOQ.get_population(_rho_ramsey_u)

	_loss = sum(abs2,_population_ramsey-population_synthetic)/N_dark_times
	return _loss
end
```
#### Step 2b: define a callback function to plot while optimizing
```julia
plot_callback = function(p,other_args)
	rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
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
```
#### Step 2c: Initial guess and bounds for the optimization
```julia
p_true = [freqs;gamma2] # values to generate synthetic data
# initial guess for the optimization
p_initial = [freqs .- 0.5e-4;0.9.*gamma2]

# bounds for the optimization
lower_bound = (0.5).*p_true
upper_bound = (1.5).*p_true
```
### Step 3: solve the optimization problem
#### Step 3a: define the optimization object (objective function and its gradient)
```julia
# construct optimization object, use Zygote auto-differentiation to compute the gradient
loss_gradient = GalacticOptim.OptimizationFunction(loss, GalacticOptim.AutoZygote())
opt_prob = GalacticOptim.OptimizationProblem(loss_gradient, p_initial,
					lb = lower_bound, ub = upper_bound)
```
#### Step 3b: solve the optimization problem with the Optim interface of `GalacticOptim.jl`
```julia
println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = GalacticOptim.solve(opt_prob ,Fminbox(LBFGS()),
			cb = plot_callback,
			outer_iterations = 20,
			iterations = 10,
			show_trace=true,
			f_tol = 1e-10,
			outer_f_tol = 1e-10)
println("Optim Fminbox(LBFGS) Optimization done")
```
#### Step 3c: solve the optimization problem but with the NLopt interface of `GalacticOptim.jl`
println("NLopt LBFGS Optimization starts")
@time sol_nlopt_LBFGS = GalacticOptim.solve(opt_prob, Opt(:LD_LBFGS,length(p_initial)),
			maxiters=200,
			cb = plot_callback,
			ftol_rel=1e-7)
println("NLopt LBFGS Optimization done")
#### Step 3d: presnet the result
```julia
# present the solutions
println("\nOptimized results: ",sol.u,
        "\nLoss: ",sol.minimum,
		"\nError: ",sol.u-p_true)
```
#### Compare optimized results and the synthetic data.
![Example 1: Optimized results v.s. Synthetic data](Example1_result.png)
##### Population of the Ramsey 0-1 experiment for a single qubit. 
- Syn-0: synthetic data for the energy level 0; 
- Syn-1: synthetic data for the energy level 1; 
- Opt-0: optimzed result for the energy level 0; 
- Opt-1: optimzed result for the energy level 1.