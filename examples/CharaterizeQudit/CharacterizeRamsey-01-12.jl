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
#data_set_option = "short"
data_set_option = "long"
do_optimization = true

# Load Ramsey 01 data
if (data_set_option=="short")
    ramsey_01_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_0.txt"
    ramsey_01_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_1.txt"
    ramsey_01_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_2.txt"
    t_ramsey_01_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_4000_5.txt"

    ramsey_12_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_confusion_0.txt"
    ramsey_12_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_confusion_1.txt"
    ramsey_12_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_confusion_2.txt"
    t_ramsey_12_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_12_1000000.0_5_4000_1000_20220413_confusion.txt"
else
    ramsey_01_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt"
    ramsey_01_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_1.txt"
    ramsey_01_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_2.txt"
    t_ramsey_01_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_80000_20.txt"

    ramsey_12_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_confusion_0.txt"
    ramsey_12_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_confusion_1.txt"
    ramsey_12_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_confusion_2.txt"
    t_ramsey_12_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/darktime_ramsey_12_1000000.0_20_80000_1000_20220413_confusion.txt"
end
detuning_01 = 1e-3

pop_R01_0 = readdlm(ramsey_01_0_file)
pop_R01_1 = readdlm(ramsey_01_1_file)
pop_R01_2 = readdlm(ramsey_01_2_file)
t_R01  = readdlm(t_ramsey_01_file)

if (do_optimization)
    N_used = min(250,length(pop_R01_0))
else
    N_used = length(pop_R01_0)
end

pop_R01_0 = pop_R01_0[2:N_used]
pop_R01_1 = pop_R01_1[2:N_used]
pop_R01_2 = pop_R01_2[2:N_used]
t_R01 = t_R01[2:N_used]

pop_R01_data = [pop_R01_0 pop_R01_1 pop_R01_2]

# Load Ramsey 12 files
#ramsey_12_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_0.txt"
#ramsey_12_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_1.txt"
#ramsey_12_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_4000_5_1000_2.txt"

#ramsey_12_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_0.txt"
#ramsey_12_1_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_1.txt"
#ramsey_12_2_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_12_1000000.0_80000_20_1000_2.txt"
detuning_12 = 1e-3

pop_R12_0 = readdlm(ramsey_12_0_file)
pop_R12_1 = readdlm(ramsey_12_1_file)
pop_R12_2 = readdlm(ramsey_12_2_file)
t_R12  = readdlm(t_ramsey_12_file)

pop_R12_0 = pop_R12_0[1:N_used]
pop_R12_1 = pop_R12_1[1:N_used]
pop_R12_2 = pop_R12_2[1:N_used]
t_R12 = t_R12[1:N_used]

pop_R12_data = [pop_R12_0 pop_R12_1 pop_R12_2]

# Set up system parameters for the open quantum system and each experiments
N_states = 3; # number of states
freqs = [3.445495118; 3.237426519] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [1.0/(160.0*GLOQ.GLOQ_MICRO_SEC); 1.0/(275.0*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
T2 = [7.5;2.5]
#T2 = [7.5;15.0]
#T2 = [200.0;200.0]
gamma2   = [1.0/(T2[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T2[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)

charge_noise12 = 1.25e-4
charge_noise01 = 0.0
charge_noise = 2.0*pi*[charge_noise01;charge_noise12]

# Ramsey 01 experiments
omr_R01 = 2.0*pi*(freqs[1]-detuning_01) # Driving frequency
TC01 = 152.0# total control time
amp_01 = 1.0

# Ramsey 12 experiments
omr_R12 = 2.0*pi*(freqs[2]-detuning_12) # Driving frequency
TC12 = 152.0# total control time
amp_12 = 1.0

# Initial state
initial_state_01 = 0
state_u0 = [0.0;0.0;0.0]
state_v0 = [0.0;0.0;0.0]
state_u0[initial_state_01+1] = 1.0

###########################################
# Results corresponding to initial guess
##########################################
# Ramsey-01
rho_syn_R01_u,_ = GLOQ.RamseyForwardSolve(
    state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
    omegas,omr_R01, # transition frequencies, drive frequency
    gamma1,gamma2, # decay and dephasing
    0, # initial state
    TC01,t_R01,N_states,
    amp_01,
    half_pi_option=half_pi_str) # control time, dark time, total number of states
pop_R01_syn = GLOQ.get_population(rho_syn_R01_u)

# Ramsey-12
rho_syn_R12_u,_ = GLOQ.PerformForwardRamseyExperimentParity(
                          state_u0,state_v0,
				   		  omegas,[omegas[1];omr_R12],
                          charge_noise,
				   		  gamma1,gamma2,
						  1,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)
pop_R12_syn = GLOQ.get_population(rho_syn_R12_u)

# Plot for Ramsey-01
fig_r01=plot(t_R01/GLOQ.GLOQ_MICRO_SEC,pop_R01_data,
             line=(2.5,:solid),
             lab = ['0' '1' '2'],
             xlabel="μs",
             ylabel="Population",
             title=string("Ramsey 0-1, detuning = ",detuning_01*1e3, "MHz") )
if (do_optimization)
    plot!(fig_r01,t_R01/GLOQ.GLOQ_MICRO_SEC,pop_R01_syn,legend=:false,
          line=(2.5,:dash))
end
display(fig_r01)

# Plot for Ramsey-12
fig_r12=plot(t_R12/GLOQ.GLOQ_MICRO_SEC,pop_R12_data,
             line=(2.5,:solid),
             lab = ['0' '1' '2'],
             xlabel="μs",
             ylabel="Population",
             title=string("Ramsey 1-2, detuning = ",detuning_12*1e3, "MHz") )
if (do_optimization)
    plot!(fig_r12,t_R12/GLOQ.GLOQ_MICRO_SEC,pop_R12_syn,legend=:false,
          line=(2.5,:dash))
end
display(fig_r12)

println("Ramsey 1-2 difference: ",sum(abs2,pop_R12_data-pop_R12_syn))

if (do_optimization)
##########
# Define the loss function for the GalacticOptim
# p: phyiscal parameters:
#    p[1,2] = transition frequencies in GHz
#    p[3,4] = gamma2s
# dummy_parameter: needed by GalacticOptim, one can just put [] here
function loss_ramsey(p,dummy_parameter)
    # Transition frequencies
    _freqs = [p[1];p[2]]
    # gamma2
    _gamma2 = [1.0/(p[3]*GLOQ.GLOQ_MICRO_SEC); 1.0/(p[4]*GLOQ.GLOQ_MICRO_SEC)]
    # charge noise
    _charge_noise = 2.0*pi*[charge_noise01;p[5]]
    # amplitude of pulse
    _amp_01 = p[6]
    _amp_12 = p[7]
    # Ramsey 0-1
    _rho_R01_u,_ = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                           (2*pi).*_freqs,omr_R01,
                                           gamma1,_gamma2,# gamma2
                                           0, # initial state
                                           TC01,t_R01,N_states,_amp_01,
                                           half_pi_option=half_pi_str)
    _pop_R01 = GLOQ.get_population(_rho_R01_u)

    # Ramsey 1-2
    _rho_R12_u,_ = GLOQ.PerformForwardRamseyExperimentParity(
                          state_u0,state_v0,
                          (2*pi).*_freqs,[omegas[1];omr_R12],
                          _charge_noise,
				   		  gamma1,_gamma2,
						  1,
				   		  [TC01;TC12],
						  [_amp_01;_amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)
    _pop_R12 = GLOQ.get_population(_rho_R12_u)

    _loss = sum(abs2,_pop_R01-pop_R01_data)+
            sum(abs2,_pop_R12-pop_R12_data)
    return _loss
end




plot_callback = function(p,other_args)
    # Transition frequencies
    _freqs = [p[1];p[2]]
    # gamma2
    _gamma2 = [1.0/(p[3]*GLOQ.GLOQ_MICRO_SEC); 1.0/(p[4]*GLOQ.GLOQ_MICRO_SEC)]
    # charge noise
    _charge_noise = 2.0*pi*[charge_noise01;p[5]]
    # amplitude of pulse
    _amp_01 = p[6]
    _amp_12 = p[7]
    # Ramsey-01
    # Ramsey 0-1
    _rho_R01_u,_ = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                           (2*pi).*_freqs,omr_R01,
                                           gamma1,_gamma2,# gamma2
                                           0, # initial state
                                           TC01,t_R01,N_states,_amp_01,
                                           half_pi_option=half_pi_str)
    _pop_R01 = GLOQ.get_population(_rho_R01_u)

    _rho_R12_u,_ = GLOQ.PerformForwardRamseyExperimentParity(
        state_u0,state_v0,
        (2*pi).*_freqs,[omegas[1];omr_R12],
        _charge_noise,
        gamma1,_gamma2,
        1,
        [TC01;TC12],
        [_amp_01;_amp_12],
        t_R12,
        half_pi_option=half_pi_str)
    _pop_R12 = GLOQ.get_population(_rho_R12_u)

    _fig_R01 = plot(t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,label=["Data-0" "Data-1" "Data-2"],
                      line = (:dash,0.0), marker = ([:hex :hex], 5, 0.1),legend=:outerright,
                      title="Ramsey 0-1");
    plot!(_fig_R01,t_R01./GLOQ.GLOQ_MICRO_SEC,_pop_R01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,line = (:solid,3.0))
    
    _fig_R12 = plot(t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,label=["Data-0" "Data-1" "Data-2"],
                      line = (:dash,0.0), marker = ([:hex :hex], 5, 0.1),legend=:outerright,
                      title="Ramsey 1-2");
    plot!(_fig_R12,t_R12./GLOQ.GLOQ_MICRO_SEC,_pop_R12,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,line = (:solid,3.0))
    ######################################
    display( plot(_fig_R01,_fig_R12,
                  legendfontsize=10,xtickfontsize=10,ytickfontsize=10,titlefontsize=12) )
    return false
end

# initial guess for the optimization
p_initial = [ freqs[1];freqs[2];T2[1];T2[2];charge_noise12;1.0;1.0 ]
# bounds for the optimization
lower_bound = [freqs[1]-1e-3;freqs[2]-1e-3;0.1; 0.1; charge_noise12-0.5e-4;0.99;0.99]#TC01*0.1]
upper_bound = [freqs[1]+1e-3;freqs[2]+1e-3;30.0;20.0;charge_noise12+0.5e-4;1.01;1.01]#TC01*5.0]

#loss_ramsey_test(x) = loss_ramsey(x,[])
#println(loss_ramsey(p_initial,[]))
#@time Zygote.gradient(loss_ramsey_test,p_initial)
#@time plot_callback(p_initial,[])

# construct optimization object, use Zygote auto-differentiation to compute the gradient
loss_gradient = GalacticOptim.OptimizationFunction(loss_ramsey, GalacticOptim.AutoZygote())
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
my_freqs  = [sol.u[1];sol.u[2]]
my_gamma2 = [1.0/(sol.u[3]*GLOQ.GLOQ_MICRO_SEC);1.0/(sol.u[4]*GLOQ.GLOQ_MICRO_SEC)]
my_charge_noise12 = 2.0*pi*sol.u[5]
my_amp_01 = sol.u[6]
my_amp_12 = sol.u[7]

# Ramsey 01 results
rho_R01_u,_ = GLOQ.RamseyForwardSolve(state_u0,state_v0,
                                      (2*pi).*my_freqs,omr_R01,
                                      gamma1,my_gamma2,# gamma2
                                      0, # initial state
                                      TC01,t_R01,N_states,my_amp_01,
                                      half_pi_option=half_pi_str)
pop_R01 = GLOQ.get_population(rho_R01_u)
fig_R01 = plot(t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title=string("Ramsey 0-1, detuning = ",detuning_01*1e3, "MHz") );
plot!(fig_R01,t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,
       line = (:dash,3.0))
display(fig_R01)


fig_R01_error = plot(t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data-pop_R01,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title=string("Ramsey 0-1, detuning = ",detuning_01*1e3, "MHz") );
display(fig_R01_error)

# Ramsey 12 results
rho_R12_u,_ = GLOQ.PerformForwardRamseyExperimentParity(
                    state_u0,state_v0,
                    (2*pi).*my_freqs,[omegas[1];omr_R12],
                    [charge_noise01;my_charge_noise12],
                    gamma1,my_gamma2,
                    1,
                    [TC01;TC12],
                    [my_amp_01;my_amp_12],
                    t_R12,
                    half_pi_option=half_pi_str)
pop_R12 = GLOQ.get_population(rho_R12_u)
fig_R12 = plot(t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title=string("Ramsey 1-2, detuning = ",detuning_12*1e3, "MHz") );
plot!(fig_R12,t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12,label=["Opt-0" "Opt-1" "Opt-2"],legend=:outerright,
       line = (:dash,3.0))
display(fig_R12)


fig_R12_error = plot(t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data-pop_R12,label=["Data-0" "Data-1" "Data-2"],
                  line = (:solid,3.0), #marker = ([:hex :hex], 5),
                  legend=:outerright,
                  xlabel="μs",
                  title=string("Ramsey 1-2 error, detuning = ",detuning_12*1e3, "MHz") );
display(fig_R12_error)

# 
println("Optimized values: ",sol.u,"\n",
        "Intitial guess:   ",p_initial,"\n",
        "Difference:       ",sol.u-p_initial)
end