using Zygote
using LinearAlgebra
using GalacticOptim,NLopt,Optim
using Plots
using GLOQ
using DelimitedFiles

pop_R01_0 = readdlm("/Users/appelo/Dropbox/Research/QuantumZhichaoDaniel/CharaterizationAndControl/data0226_ramsey_01/ramsey_01_250000.0_1000_0.txt")
pop_R01_1 = readdlm("/Users/appelo/Dropbox/Research/QuantumZhichaoDaniel/CharaterizationAndControl/data0226_ramsey_01/ramsey_01_250000.0_1000_1.txt")
pop_R01_2 = readdlm("/Users/appelo/Dropbox/Research/QuantumZhichaoDaniel/CharaterizationAndControl/data0226_ramsey_01/ramsey_01_250000.0_1000_2.txt")

# System parameters for a simple two level open quantum system
N_states = 3; # number of states
freqs = [3.4511896626; 3.242837221] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [1.0/(45.0*GLOQ.GLOQ_MICRO_SEC); 1.0/(45.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
gamma2   = [1.0/(24.0*GLOQ.GLOQ_MICRO_SEC); 1.0/(24.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)
omr_R01 = 2.0*pi*(3.4512896626 - 2.5e-4) # drive frequency for the Ramsey experiment
TC01 = 2.5*100.0 # total control time
# Initial state
initial_state = 0
state_u0 = [0.0;0.0;0.0]
state_v0 = [0.0;0.0;0.0]
state_u0[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time
T_R01 = 1000*20
# total number of dark time samples
N_dark_R01 = 1001

dt_R01 = T_R01/N_dark_R01
t_dark_R01 = collect(range(0.0, T_R01, length=N_dark_R01))
# Forward solve
rho_syn_R01_u,rho_syn_R01_v = GLOQ.RamseyForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas,omr_R01, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    0, # initial state
    TC01,t_dark_R01,N_states) # control time, dark time, total number of states

pop_R01_syn = GLOQ.get_population(rho_syn_R01_u)

fig01=plot(t_dark_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_syn)

omr_R12 = 2.0*pi*(3.242837221 - 5.0e-4) # drive frequency for the Ramsey experimen
TC12 = 2.5*100.0 # total control time
# Initial state
initial_state = 1
state_u1 = [0.0;0.0;0.0]
state_v1 = [0.0;0.0;0.0]
state_u1[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time
T_R12 = 1000*20
# total number of dark time samples
N_dark_R12 = 1001

dt_R12 = T_R12/N_dark_R12
t_dark_R12 = collect(range(0.0, T_R12, length=N_dark_R12))
# Forward solve
rho_syn_R12_u,rho_syn_R12_v = GLOQ.RamseyForwardSolve(
    state_u1,state_v1, # initial values, u for the real part, v for the imaginary part
    omegas,omr_R12, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    1, # initial state
    TC12,t_dark_R12,N_states) # control time, dark time, total number of states

pop_R12_syn = GLOQ.get_population(rho_syn_R12_u)

fig12=plot(t_dark_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_syn)

# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#    p[1,2] = transition frequencies in GHz
#    p[3,4] = gamma2s
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss(p,dummy_parameter)
    # Ramsey
    _rho_R01_u,_rho_R01_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                                    (2*pi).*[p[1];p[2]],omr_R01,
                                                    gamma1,[p[3];p[4]],# gamma2
                                                    0, # initial state
                                                    TC01,t_dark_R01,N_states)
    _pop_R01 = GLOQ.get_population(_rho_R01_u)

    _loss = sum(abs2,_pop_R01-pop_R01_syn)*dt_R01#+
#        sum(abs2,_population_echo-population_echo_synthetic)*dt_echo+
#        sum(abs2,_population_t1-population_t1_synthetic)*dt_t1

    return _loss
end


plot_callback = function(p,other_args)
    rho_R01_u,rho_R01_v = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                                  (2*pi).*[p[1];p[2]],omr_R01,
                                                  gamma1,[p[3];p[4]],# gamma2
                                                  0, # initial state
                                                  TC01,t_dark_R01,N_states)
    pop_R01 = GLOQ.get_population(rho_R01_u)
    # Ramsey
    fig_R01 = plot(t_dark_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_syn,label=["Syn-0" "Syn-1" "Syn-2"],
                      line = (:dash,0.0), marker = ([:hex :hex], 5, 0.5),legend=:outerright,
                      title="Ramsey");
    plot!(fig_R01,t_dark_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright);
    display(fig_R01)
    #    display( plot(fig_ramsey,fig_echo,fig_t1,layout=grid(3,1),size=[1000,1500],
#                  legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18) )
    return false
end
=#

# initial guess for the optimization
p_initial = [freqs.-0.5e-4;1.25.*gamma2]
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
                                f_tol = 1e-5,
                                outer_f_tol = 1e-5)
println("Optim Fminbox(LBFGS) Optimization done")

# present the solutions
println("\nOptimized results: ",sol.u,
        "\nLoss: ",sol.minimum,
        "\nError: ",sol.u-p_true)
=#
