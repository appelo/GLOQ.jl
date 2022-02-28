using GLOQ, Zygote, LinearAlgebra
using GalacticOptim,NLopt,Optim
using Plots
using DelimitedFiles

pop_0 = readdlm("data0226_ramsey/ramsey_01_250000.0_1000_0.txt")
pop_1 = readdlm("data0226_ramsey/ramsey_01_250000.0_1000_1.txt")

# System parameters for a simple two level open quantum system
N_states = 2; # number of states
freqs = [3.4511896626] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [1.0/(55.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
gamma2   = [1.0/(15.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)
omr = 2.0*pi*(3.4512896626-2.5e-4) # drive frequency
TC = 2.5*100.0 # total control time in ns needed for a Rabi pulse

# Initial state
initial_state = 0
state_u0 = [0.0;0.0]
state_v0 = [0.0;0.0]
state_u0[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time given in Microseconds
T_Ramsey = 1000*20 # ns #40.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_times = 1001
t_dark_times = collect(range(0.0, T_Ramsey, length=N_dark_times))

# Forward solve to generate synthetic data
rho_synthetic_ramsey_u,rho_synthetic_ramsey_v = GLOQ.RamseyForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas,omr, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing parameters
    initial_state, # initial state
    TC,t_dark_times,N_states) # control time, dark time, total number of states
population_synthetic = GLOQ.get_population(rho_synthetic_ramsey_u)

pop_data = [pop_0 pop_1]

pl1 = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
plot!(pl1,t_dark_times./GLOQ.GLOQ_MICRO_SEC,pop_data)

display(pl1)

population_synthetic = pop_data


# Plot the synthetic data
fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
display(fig)

# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#        p[1] = transition frequency in GHz
#        p[2] = gamma2
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p,dummy_parameter)
    _rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                                          (2*pi).*[p[1]],omr,
                                                          gamma1,[p[2]],#gamma1,gamma2,
                                                          initial_state, # initial state
                                                          TC,t_dark_times,N_states)
    _population_ramsey = GLOQ.get_population(_rho_ramsey_u)

    _loss = sum(abs2,_population_ramsey-population_synthetic)/N_dark_times
    return _loss
end

plot_callback = function(p,other_args)
    rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                                        (2*pi).*[p[1]],omr,
                                                        gamma1,[p[2]],#gamma1,gamma2,
                                                        initial_state, # initial state
                                                        TC,t_dark_times,N_states)
    population_ramsey = GLOQ.get_population(rho_ramsey_u)
    fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic,label=["Syn-0" "Syn-1"],
             line = (:dash,0.0), marker = ([:hex :hex], 5, 0.05)  )
    plot!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_ramsey,label=["Opt-0" "Opt-1"]
          ,line = (:solid, 4))                     #,line = (:dot, 4), marker = ([:hex :hex], 5, 0.1))
    display(fig)
    return false
end


p_true = [3.4511896626;gamma2] # values to generate synthetic data
# initial guess for the optimization
p_initial = [3.4511896626;0.9.*gamma2]

# bounds for the optimization
lower_bound = [3.4510;gamma2*0.1]
upper_bound = [3.4513;gamma2*2.5]


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
                                f_tol = 1e-5,
                                outer_f_tol = 1e-5)
println("Optim Fminbox(LBFGS) Optimization done")
# present the solutions
println("\nOptimized results: ",sol.u,
        "\nLoss: ",sol.minimum)
#=
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
=#
