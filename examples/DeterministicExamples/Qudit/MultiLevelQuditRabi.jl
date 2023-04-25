using Zygote, LinearAlgebra
using Optim
using Plots
using DelimitedFiles
using GLOQ


# Read the experimental data
pop_R01_0 = readdlm("population20220321/population_rabi_01_1000_2432_0_20220321_0.txt")
pop_R01_1 = readdlm("population20220321/population_rabi_01_1000_2432_0_20220321_1.txt")
pop_R01_2 = readdlm("population20220321/population_rabi_01_1000_2432_0_20220321_2.txt")
# Stack the data
pop_R01_data = [pop_R01_0 pop_R01_1 pop_R01_2]

# Array containg all the times when the populations are recorded.
t_rabi = readdlm("population20220321/time_rabi_01_1000_2432_0_20220321.txt")

# System parameters for a simple two level open quantum system
N_states = 3; # number of states
freqs = [3.451266428-1e-4; 3.242951102] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
# The T1-T2 times should not matter much for this example. We choose them to be something that is reasonable. 
gamma1   = [1.0/(45.0*GLOQ.GLOQ_MICRO_SEC); 1.0/(275.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
T2 = [15.0;20.0]
gamma2   = [1.0/(T2[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T2[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)
omr_01 = 2.0*pi*(3.451266428) #freqs[1]) # drive frequency
TC01 = 152.0 # total control time in ns needed for a Rabi pulse 
amp_01 = 1.0

# Initial state
initial_state_01 = 0
state_u0 = [0.0;0.0;0.0]
state_v0 = [0.0;0.0;0.0]
state_u0[initial_state_01+1] = 1.0

##########
# Forward solve
rho_syn_R01_u,_ = GLOQ.RabiForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas,omr_01, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    initial_state_01, # initial state
    TC01,t_rabi,N_states,amp_01) # control time, dark time, total number of states
pop_R01_syn = GLOQ.get_population(rho_syn_R01_u)

rho_syn_R01_u2,_ = GLOQ.RabiForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas.-2*pi*1e-4,omr_01, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    initial_state_01, # initial state
    TC01,t_rabi,N_states,amp_01) # control time, dark time, total number of states
pop_R01_syn2 = GLOQ.get_population(rho_syn_R01_u2)

fig_r01=plot(t_rabi,pop_R01_data,
             line=(0.0,:solid),
             marker=(:hex))

plot!(fig_r01,t_rabi,(sin.(0.5*pi/154*t_rabi)).^2,legend=:false,
      line=(2.5,:solid), 
      xlabel="ns",
      ylabel="Population")


fig=plot(t_rabi,pop_R01_syn,legend=:false,
      line=(2.5,:solid), 
      xlabel="ns",
      ylabel="Population")
plot!(t_rabi,pop_R01_syn2,legend=:false,
      line=(2.5,:dash), 
      xlabel="ns",
      ylabel="Population")
display(fig)



# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#	 p[1] = transition frequency in GHz
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p)
    _amp = p[2]
    _rho_R01_u,_rho_R01_v = GLOQ.RabiForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        [(2*pi)*p[1]; omegas[2]],omr_01, # transition frequencies, drive frequency
        gamma1,gamma2, # decay and dephasing
        initial_state_01,
        TC01,t_rabi,N_states,_amp) # control time, dark time, total number of states
    _pop_R01 = GLOQ.get_population(_rho_R01_u)
    _loss = sum(abs2,_pop_R01-pop_R01_data)
    return _loss
end

function loss_gradient!(G,p)
    G .= Zygote.gradient(loss,p)[1]
end

function loss_nlopt(p,grad)
    grad .= Zygote.gradient(loss,p)[1]
    return loss(p)
end

plot_callback = function(_trace)
    p = _trace.metadata["x"]
    _amp = p[2]
    rho_R01_u,_rho_R01_v = GLOQ.RabiForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        [(2*pi)*p[1]; omegas[2]],omr_01, # transition frequencies, drive frequency
        gamma1,gamma2, # decay and dephasing
        initial_state_01,
        TC01,t_rabi,N_states,_amp) # control time, dark time, total number of states
    pop_R01 = GLOQ.get_population(rho_R01_u)
    fig=plot(t_rabi,pop_R01_data,label=["Data-0" "Data-1" "Data-2"],
	     line = (:dash,0.0), marker = ([:hex :hex :hex], 5, 0.5),
         xlabel = "Duration of the pulse (ns)",ylabel = "Population"  )
    plot!(fig,t_rabi,pop_R01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright
	  )
    display(fig)
    return false
end

# initial guess for the optimization
p_initial = [freqs[1]*1.0001;amp_01]

# bounds for the optimization
lower_bound = [p_initial[1]-5e-3;0.5]
upper_bound = [p_initial[1]+5e-3;1.5]

println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = Optim.optimize(loss,loss_gradient!,
                lower_bound,upper_bound,p_initial,
                Fminbox(LBFGS()),
                Optim.Options(
				callback = plot_callback,extended_trace=true,
				outer_iterations = 20,
				iterations = 10,
				show_trace=true,
				f_tol = 1e-10,
				outer_f_tol = 1e-10))
println("Optim Fminbox(LBFGS) Optimization done")
println("\nOptim results: ",sol.minimizer,
        "\nOptim loss: ",sol.minimum)

println("NLopt LBFGS Optimization starts")
nlopt_obj = NLopt.Opt(:LD_LBFGS,length(p_initial))
nlopt_obj.min_objective = loss_nlopt
nlopt_obj.lower_bounds = lower_bound
nlopt_obj.upper_bounds = upper_bound
nlopt_obj.maxeval = 200
nlopt_obj.ftol_rel = 1e-7
@time loss_value_nlopt,sol_nlopt,flag_nlopt = NLopt.optimize(nlopt_obj,p_initial)
println("NLopt LBFGS Optimization done")

println("\nNLopt results: ",sol_nlopt,
        "\nNLopt loss: ",loss_value_nlopt)






