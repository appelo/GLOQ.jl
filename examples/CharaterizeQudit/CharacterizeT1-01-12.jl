include("../../src/GLOQ.jl")
using Zygote
using LinearAlgebra
using GalacticOptim,NLopt,Optim
using Plots
#using GLOQ
using DelimitedFiles

half_pi_str = "duration" 
#############################################################################
# Read data files
#############################################################################
#t1_01_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_t1_01_1000_20220413_40000_0.txt"
t1_01_0_file = "data-set-20220413/population_t1_01_1000_20220413_40000_0.txt"
t1_01_1_file = "data-set-20220413/population_t1_01_1000_20220413_40000_1.txt"
t1_01_2_file = "data-set-20220413/population_t1_01_1000_20220413_40000_2.txt"
t_t1_01_file = "data-set-20220413/darktime_t1_01_1000_20220413_40000.txt"

t1_12_0_file = "data-set-20220413/population_t1_12_1000_40000_20220413_0.txt"
t1_12_1_file = "data-set-20220413/population_t1_12_1000_40000_20220413_1.txt"
t1_12_2_file = "data-set-20220413/population_t1_12_1000_40000_20220413_2.txt"
t_t1_12_file = "data-set-20220413/darktime_t1_12_1000_40000_20220413.txt"

pop_T01_0 = readdlm(t1_01_0_file)
pop_T01_1 = readdlm(t1_01_1_file)
pop_T01_2 = readdlm(t1_01_2_file)
t_T01  = readdlm(t_t1_01_file)

pop_T01_0 = pop_T01_0[8:end]
pop_T01_1 = pop_T01_1[8:end]
pop_T01_2 = pop_T01_2[8:end]
t_T01 = t_T01[8:end]

pop_T01_data = [pop_T01_0 pop_T01_1 pop_T01_2]

# Load T1 12 files
pop_T12_0 = readdlm(t1_12_0_file)
pop_T12_1 = readdlm(t1_12_1_file)
pop_T12_2 = readdlm(t1_12_2_file)
t_T12  = readdlm(t_t1_12_file)

pop_T12_data = [pop_T12_0 pop_T12_1 pop_T12_2]

# Set up system parameters for the open quantum system and each experiments
N_states = 3; # number of states
freqs = [3.445495118; 3.237426519] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
T1 = [160.0;100.0]
gamma1   = [1.0/(T1[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T1[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
T2 = [7.5;2.5]
#T2 = [7.5;15.0]
#T2 = [200.0;200.0]
gamma2   = [1.0/(T2[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T2[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)

charge_noise12 = 1.25e-4
charge_noise01 = 0.0
charge_noise = 2.0*pi*[charge_noise01;charge_noise12]

# Ramsey 01 experiments
omr_T01 = 2.0*pi*(freqs[1]) # Driving frequency
TC01 = 152.0# total control time
amp_01 = 1.0

# Ramsey 12 experiments
omr_T12 = 2.0*pi*(freqs[2]) # Driving frequency
TC12 = 152.0# total control time
amp_12 = 1.0

# Initial state
state_u0 = [1.0;0.0;0.0]
state_v0 = [0.0;0.0;0.0]

state_u1 = [0.0;1.0;0.0]
state_v1 = [0.0;0.0;0.0]
###########################################
# Results corresponding to initial guess
##########################################
# T1-01
rho_syn_R01_u,_ = GLOQ.T1ForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas,omr_T01, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    0, # initial state
    TC01,t_T01,N_states) # control time, dark time, total number of states
pop_T01_syn = GLOQ.get_population(rho_syn_R01_u)

# T1-12
rho_syn_T12_u,_ = GLOQ.T1ForwardSolve(
    state_u1,state_v1, # initial values, u for the real part, v for the imaginary part
    omegas,omr_T12, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    1, # initial state
    TC01,t_T12,N_states)
pop_T12_syn = GLOQ.get_population(rho_syn_T12_u)

# Plot for Ramsey-01
fig_t01=plot(t_T01/GLOQ.GLOQ_MICRO_SEC,pop_T01_data,
             line=(2.5,:solid),
             lab = ['0' '1' '2'],
             xlabel="μs",
             ylabel="Population",
             title="T1 0-1" )
plot!(fig_t01,t_T01/GLOQ.GLOQ_MICRO_SEC,pop_T01_syn,legend=:false,
      line=(2.5,:dash))
display(fig_t01)

# Plot for Ramsey-12
fig_t12=plot(t_T12/GLOQ.GLOQ_MICRO_SEC,pop_T12_data,
             line=(2.5,:solid),
             lab = ['0' '1' '2'],
             xlabel="μs",
             ylabel="Population",
             title="T1 1-2" )
plot!(fig_t12,t_T12/GLOQ.GLOQ_MICRO_SEC,pop_T12_syn,legend=:false,
      line=(2.5,:dash))
display(fig_t12)

##########
# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#    p[1,2] = transition frequencies in GHz
#    p[3,4] = gamma2s
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss_T1(p,dummy_parameter)
    # gamma1
    _gamma1 = [1.0/(p[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(p[2]*GLOQ.GLOQ_MICRO_SEC)]

    # T1 0-1
    _rho_T01_u,_ = GLOQ.T1ForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T01, # transition frequencies, drive frequency
        _gamma1,gamma2, # decay and dephasing
        0, # initial state
        TC01,t_T01,N_states) # control time, dark time, total number of states
    _pop_T01 = GLOQ.get_population(_rho_T01_u)

    # T1 1-2
    _rho_T12_u,_ = GLOQ.T1ForwardSolve(
        state_u1,state_v1, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T12, # transition frequencies, drive frequency
        _gamma1,gamma2, # decay and dephasing
        1, # initial state
        TC12,t_T12,N_states) # control time, dark time, total number of states
    _pop_T12 = GLOQ.get_population(_rho_T12_u)

    _loss = sum(abs2,_pop_T01-pop_T01_data)+
            sum(abs2,_pop_T12-pop_T12_data)
    return _loss
end


plot_callback = function(p,other_args)
    # gamma1
    _gamma1 = [1.0/(p[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(p[2]*GLOQ.GLOQ_MICRO_SEC)]

    # T1 0-1
    _rho_T01_u,_ = GLOQ.T1ForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T01, # transition frequencies, drive frequency
        _gamma1,gamma2, # decay and dephasing
        0, # initial state
        TC01,t_T01,N_states) # control time, dark time, total number of states
    _pop_T01 = GLOQ.get_population(_rho_T01_u)

    # T1 1-2
    _rho_T12_u,_ = GLOQ.T1ForwardSolve(
        state_u1,state_v1, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T12, # transition frequencies, drive frequency
        _gamma1,gamma2, # decay and dephasing
        1, # initial state
        TC12,t_T12,N_states) # control time, dark time, total number of states
    _pop_T12 = GLOQ.get_population(_rho_T12_u)

    _fig_T01 = plot(t_T01./GLOQ.GLOQ_MICRO_SEC,pop_T01_data,label=["Data-0" "Data-1" "Data-2"],
                      line = (:dash,0.0), marker = ([:hex :hex], 5, 0.1),legend=:outerright,
                      title="T1 0-1");
    plot!(_fig_T01,t_T01./GLOQ.GLOQ_MICRO_SEC,_pop_T01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,line = (:solid,3.0))
    
    _fig_T12 = plot(t_T12./GLOQ.GLOQ_MICRO_SEC,pop_T12_data,label=["Data-0" "Data-1" "Data-2"],
                      line = (:dash,0.0), marker = ([:hex :hex], 5, 0.1),legend=:outerright,
                      title="T1 1-2");
    plot!(_fig_T12,t_T12./GLOQ.GLOQ_MICRO_SEC,_pop_T12,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,line = (:solid,3.0))
    ######################################
    display( plot(_fig_T01,_fig_T12,
                  legendfontsize=10,xtickfontsize=10,ytickfontsize=10,titlefontsize=12) )
    return false
end

# initial guess for the optimization
p_initial = [ T1[1];T1[2] ]
# bounds for the optimization
lower_bound = [50;75.0]#TC01*0.1]
upper_bound = [200.0;150.0]#TC01*5.0]

loss_T1_test(x) = loss_T1(x,[])
println(loss_T1(p_initial,[]))
@time Zygote.gradient(loss_T1_test,p_initial)
@time plot_callback(p_initial,[])

# construct optimization object, use Zygote auto-differentiation to compute the gradient
loss_gradient = GalacticOptim.OptimizationFunction(loss_T1, GalacticOptim.AutoZygote())
opt_prob = GalacticOptim.OptimizationProblem(loss_gradient, p_initial,
                                             lb = lower_bound, ub = upper_bound)


println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = GalacticOptim.solve(opt_prob ,Fminbox(LBFGS()),
                                cb = plot_callback,
                                outer_iterations = 50,
                                iterations = 50,
                                show_trace=true,
                                f_tol = 1e-4,
                                outer_f_tol = 1e-4)
println("Optim Fminbox(LBFGS) Optimization done")

# present the solutions
println("\nOptimized results: ",sol.u,
        "\nLoss: ",sol.minimum)

# present results
my_gamma1 = [1.0/(sol.u[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(sol.u[2]*GLOQ.GLOQ_MICRO_SEC)]
# T1 01 results
rho_T01_u,_ = GLOQ.T1ForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T01, # transition frequencies, drive frequency
        my_gamma1,gamma2, # decay and dephasing
        0, # initial state
        TC01,t_T01,N_states) # control time, dark time, total number of states
pop_T01 = GLOQ.get_population(rho_T01_u)

fig_T01 = plot(t_T01./GLOQ.GLOQ_MICRO_SEC,pop_T01_data,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title="T1 0-1" );
plot!(fig_T01,t_T01./GLOQ.GLOQ_MICRO_SEC,pop_T01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,
       line = (:dash,3.0))
display(fig_T01)


fig_T01_error = plot(t_T01./GLOQ.GLOQ_MICRO_SEC,pop_T01_data-pop_T01,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title="T1 0-1 error" );
display(fig_T01_error)

# T1 1-2 results
rho_T12_u,_ = GLOQ.T1ForwardSolve(
        state_u1,state_v1, # initial values, u for the real part, v for the imaginary part
        omegas,omr_T12, # transition frequencies, drive frequency
        my_gamma1,gamma2, # decay and dephasing
        1, # initial state
        TC12,t_T12,N_states) # control time, dark time, total number of states
pop_T12 = GLOQ.get_population(rho_T12_u)

fig_T12 = plot(t_T12./GLOQ.GLOQ_MICRO_SEC,pop_T12_data,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title="T1 1-2" );
plot!(fig_T12,t_T12./GLOQ.GLOQ_MICRO_SEC,pop_T12,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,
       line = (:dash,3.0))
display(fig_T12)


fig_T12_error = plot(t_T12./GLOQ.GLOQ_MICRO_SEC,pop_T12_data-pop_T12,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title="T1 1-2 error" );
display(fig_T12_error)

# 
println("Optimized values: ",sol.u,"\n",
        "Intitial guess:   ",p_initial,"\n",
        "Difference:       ",sol.u-p_initial)
