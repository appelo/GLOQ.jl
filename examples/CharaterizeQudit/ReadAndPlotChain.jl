include("../../src/GLOQ.jl")
using LinearAlgebra
using Zygote
using Turing, Distributions, DifferentialEquations
# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots
using StatsPlots
using CSV,DataFrames
# Set a seed for reproducibility.
using Random
using DelimitedFiles
Random.seed!(16);
using Plots
using Dates
using StatsBase
using Distributions
using Printf
# Read data
my_chain = read("data-2022-04-15/chain_2500.jls",Chains)
my_chain_data = DataFrame(my_chain)
my_chain_size = 2500
# Extract chains for essential values
freq_01_chain = my_chain_data._freq_01
freq_12_chain = my_chain_data._freq_12
char_12_chain = my_chain_data._char_12

freq_01_mean = mean(freq_01_chain)
freq_12_mean = mean(freq_12_chain)
char_12_mean = mean(char_12_chain)

# plot the chain
fig_chain = plot(my_chain)
range_freq_01 = maximum(freq_01_chain)-minimum(freq_01_chain)
xticks!(fig_chain[6],
        [freq_01_mean-6e-7;freq_01_mean;freq_01_mean+6e-7],
        ["mean - 0.6 KHz";string(@sprintf("%.6f",freq_01_mean)," GHz");"mean + 0.6 KHz"]);
yticks!(fig_chain[5],
        [freq_01_mean-6e-7;freq_01_mean;freq_01_mean+6e-7],
        ["mean - 0.6 KHz";string(@sprintf("%.6f",freq_01_mean)," GHz");"mean + 0.6 KHz"]);
range_freq_12 = maximum(freq_12_chain)-minimum(freq_12_chain)
xticks!(fig_chain[8],
        [freq_12_mean-1.0e-6;freq_12_mean;freq_12_mean+1.0e-6],
        ["mean - 1 KHz";string(@sprintf("%.6f",freq_12_mean)," GHz");"mean + 1 KHz"]);
yticks!(fig_chain[7],
        [freq_12_mean-1.0e-6;freq_12_mean;freq_12_mean+1.0e-6],
        ["mean - 1 KHz";string(@sprintf("%.6f",freq_12_mean)," GHz");"mean + 1 KHz"]);
range_char_12 = maximum(char_12_chain)-minimum(char_12_chain)
xticks!(fig_chain[10],
       [char_12_mean-1e-6;char_12_mean;char_12_mean+1.0e-6],
       ["mean - 1 KHz";string(@sprintf("%.6f",char_12_mean*1e3)," MHz");"mean + 1 KHz"]);
yticks!(fig_chain[9],
       [char_12_mean-1e-6;char_12_mean;char_12_mean+1.0e-6],
       ["mean - 1 KHz";string(@sprintf("%.6f",char_12_mean*1e3)," MHz");"mean + 1 KHz"]);
display(fig_chain)
#display(fig_chain[7:8])
#display(fig_chain[9:10])

freqs_opt = [3.4454884101594407;3.2374500663104677]
prior_01 = TruncatedNormal(freqs_opt[1],5e-5,freqs_opt[1]-5e-4,freqs_opt[1]+5e-4)
hist_01 = fit(Histogram,freq_01_chain,nbins=500)

xmin = freqs_opt[1]-5e-4
xmax = freqs_opt[1]+5e-4

fig_hist_01 = plot(hist_01,
                   label=:false,title="Histogram for 0-1 frequency")
plot!(x->30/8000*pdf(prior_01, x),line=(:solid,5),
      xlim=[minimum(freq_01_chain),maximum(freq_01_chain)],
      label="Rescaled Prior",
      legend=:right)
xticks!(fig_hist_01,[freq_01_mean-6e-7;freq_01_mean;freq_01_mean+6e-7],    
        ["mean - 0.6 KHz";string(@sprintf("%.6f",freq_01_mean)," GHz");"mean + 0.6 KHz"]);
display(fig_hist_01)
      #xlim=[xmin,xmax])

prior_12 = truncated(Normal(freqs_opt[2],5e-5),freqs_opt[2]-2.5e-4,freqs_opt[2]+2.5e-4)
hist_12 = fit(Histogram,freq_12_chain,nbins=500)
fig_hist_12 = plot(hist_12,
                   title="Histogram for 1-2 frequency",label=:false)
xticks!(fig_hist_12,[freq_12_mean-1.0e-6;freq_12_mean;freq_12_mean+1.0e-6],
                    ["mean - 1 KHz";string(@sprintf("%.6f",freq_12_mean)," GHz");"mean + 1 KHz"]);
                    #["mean - 1 KHz";string(freq_12_mean," GHz");"mean + 1 KHz"]);
plot!(fig_hist_12,x->25/8000*pdf(prior_12, x),
      line=(:solid,5),
      xlim=[minimum(freq_12_chain),maximum(freq_12_chain)],
      label="Rescaled Prior",
      legend=:right)
display(fig_hist_12)


hist_char_12 = fit(Histogram,char_12_chain,nbins=500)
charge_noise_12 = 1.2305427472144914e-4
prior_char_12 = truncated(Normal(charge_noise_12,0.125e-5),charge_noise_12-2.5e-5,charge_noise_12+2.5e-5)
fig_hist_char_12 = plot(hist_char_12,label=:false,
                        title="Histogram for charge noise 1-2")
xticks!(fig_hist_char_12,[char_12_mean-1e-6;char_12_mean;char_12_mean+1.0e-6],
                    ["mean - 1 KHz";string(@sprintf("%.6f",char_12_mean*1e3)," MHz");"mean + 1 KHz"]);
plot!(fig_hist_char_12,x->20.0/3e5*pdf(prior_char_12, x),
      line=(:solid,5),
      xlim=[minimum(char_12_chain),maximum(char_12_chain)],
      label="Rescaled Prior",
      legend=:topright)
display(fig_hist_char_12)
#############################################################################
# plot the forward solve results 
#############################################################################
# Read data files
half_pi_str = "duration" 
data_set_option = "short"
#data_set_option = "long"

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

N_used = min(250,length(pop_R01_0))

pop_R01_0 = pop_R01_0[2:N_used]
pop_R01_1 = pop_R01_1[2:N_used]
pop_R01_2 = pop_R01_2[2:N_used]
t_R01 = t_R01[2:N_used]

pop_R01_data = [pop_R01_0 pop_R01_1 pop_R01_2]

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
freqs = [3.445495118; 3.237426519] # transition frequency used on the device in GHz
freqs_opt = [3.4454884101594407;3.2374500663104677] # transition frequency charaterized with deterministic code
omegas = 2.0*pi.*freqs # change to angular frequency
T1 = [169.1;85.79]
gamma1   = [1.0/(T1[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T1[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
#T2 = [7.5;2.5]
T2 = [6.9309771664384865;1.965597856818199]
charge_noise_01 = 0.0
charge_noise_12 = 1.2305427472144914e-4
gamma2   = [1.0/(T2[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T2[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)

charge_noise = 2.0*pi*[charge_noise_01;charge_noise_12]

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

rho_chain_mean_01,_ = GLOQ.RamseyForwardSolve(
				 state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
				 2.0*pi*[freq_01_mean;freq_12_mean],omr_R01, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing parameters ?
				 initial_state_01, # initial state
				 TC01,t_R01,N_states;
				 method="exponential")
population_chain_mean_01 = GLOQ.get_population(rho_chain_mean_01)
fig_mean_vs_noisy_01 = plot(t_R01./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_01,
		   		  label=["Mean-0" "Mean-1" "Mean-2"],
                  line=(:solid,2.5))
scatter!(fig_mean_vs_noisy_01,t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,
		 label=["Data-0" "Data-1" "Data-2"],legend=:outerright,
		 title = "Ramsey 0-1",
		 xlabel= "μs")
display(fig_mean_vs_noisy_01)

rho_chain_mean_12,_ = GLOQ.PerformForwardRamseyExperimentParity(
                          state_u0,state_v0,
                          2.0*pi*[freq_01_mean;freq_12_mean],[omegas[1];omr_R12],
                          (2.0*pi).*[charge_noise_01;char_12_mean],
				   		  gamma1,gamma2,
						  1,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)
population_chain_mean_12 = GLOQ.get_population(rho_chain_mean_12)
fig_mean_vs_noisy_12 = plot(t_R12./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_12,
		   		  line=(:solid,2.5),label=["Mean-0" "Mean-1" "Mean-2"])
scatter!(fig_mean_vs_noisy_12,t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,
		 label=["Data-0" "Data-1" "Data-2"],legend=:outerright,
		 title = "Ramsey 1-2",
		 xlabel= "μs")
display(fig_mean_vs_noisy_12)


#########################
# Sample from the chain
#########################

global fig_result_01
global fig_result_12
num_random_sample = 200
for i = 1:num_random_sample
	global fig_result_01
	global fig_result_12
	sample_ind = rand(1:my_chain_size)
    freq_01_sample = my_chain_data._freq_01[sample_ind]
    freq_12_sample = my_chain_data._freq_12[sample_ind]
    char_12_sample = my_chain_data._char_12[sample_ind]

    rho_sample_u_01,_ = GLOQ.RamseyForwardSolve(
        state_u0,state_v0, # initial values, u for the real part, v for the imaginary part
        (2.0*pi).*[freq_01_sample;freq_12_sample],omr_R01, # transition frequencies, drive frequency
        gamma1,gamma2, # decay and dephasing
        0, # initial state
        TC01,t_R01,N_states,
        amp_01,
        half_pi_option=half_pi_str)
    population_sample_01 = GLOQ.get_population(rho_sample_u_01)

	rho_sample_u_12,_ = GLOQ.PerformForwardRamseyExperimentParity(
                          state_u0,state_v0,
				   		  (2.0*pi).*[freq_01_sample;freq_12_sample],
                          [omegas[1];omr_R12],
                          (2.0*pi).*[charge_noise_01;char_12_sample],
				   		  gamma1,gamma2,
						  1,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)
	population_sample_12 = GLOQ.get_population(rho_sample_u_12)

	if (i==1)
    	fig_result_01=plot(t_R01./GLOQ.GLOQ_MICRO_SEC,population_sample_01,
              color = "#BBBBBB", label="");
		fig_result_12=plot(t_R12./GLOQ.GLOQ_MICRO_SEC,population_sample_12,
              color = "#BBBBBB", label="");
	else
    	plot!(fig_result_01,t_R01./GLOQ.GLOQ_MICRO_SEC,population_sample_01,
              color = "#BBBBBB", label="");
		plot!(fig_result_12,t_R12./GLOQ.GLOQ_MICRO_SEC,population_sample_12,
              color = "#BBBBBB", label="");
	end
end
plot!(fig_result_01,t_R01./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_01,label="",
	  xlabel="μs",
	  title=string("Ramsey 0-1, detuning = ",detuning_01*1e3," MHz")
	  )
plot!(fig_result_01,t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,label="")
display(fig_result_01)

plot!(fig_result_12,t_R12./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_12,label="",
	  xlabel="μs",
	  title=string("Ramsey 1-2, detuning = ",detuning_12*1e3," MHz")
		   		  )
plot!(fig_result_12,t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,label="")
display(fig_result_12)



#chain_data_frame = CSV.read("data-2022-04-18/chain_20_5.csv", DataFrame)
# resume sampling
#@time chain = sample(model, NUTS(0.65), sample_chain_size,
#                     save_state=true,resume_from = my_chain)