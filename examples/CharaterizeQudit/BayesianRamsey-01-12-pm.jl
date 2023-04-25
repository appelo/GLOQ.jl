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

#############################################################################
# (1) Read data files 
# (2) Set up the Turing.jl backend and define the @model macro for Bayesian inference
#############################################################################
include("ReadDataAndDefineModel-BayesianRamsey-01-12-pm.jl")

# save the chain or not
save_file = false
###########################################
# Results corresponding to initial guess
##########################################
# Ramsey-01
rho_syn_R01_u,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  omegas_minus,omegas_plus,
                          [omr_R01;omegas[2]],
				   		  gamma1,gamma2,
						  0,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R01,
                          half_pi_option=half_pi_str) # control time, dark time, total number of states
pop_R01_syn = GLOQ.get_population(rho_syn_R01_u)

# Ramsey-12
rho_syn_R12_u,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  omegas_minus,omegas_plus,
                          [omegas[1];omr_R12],
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
plot!(fig_r01,t_R01/GLOQ.GLOQ_MICRO_SEC,pop_R01_syn,legend=:false,
      line=(2.5,:dash))
display(fig_r01)

# Plot for Ramsey-12
fig_r12=plot(t_R12/GLOQ.GLOQ_MICRO_SEC,pop_R12_data,
             line=(2.5,:solid),
             lab = ['0' '1' '2'],
             xlabel="μs",
             ylabel="Population",
             title=string("Ramsey 1-2, detuning = ",detuning_12*1e3, "MHz") )
plot!(fig_r12,t_R12/GLOQ.GLOQ_MICRO_SEC,pop_R12_syn,legend=:false,
      line=(2.5,:dash))
display(fig_r12)

println("Ramsey 1-2 difference: ",sum(abs2,pop_R12_data-pop_R12_syn))

###########################################
# Bayesian sampling with NUTS sampler
##########################################
model = RamseyExperiment_01_12([pop_R01_data;pop_R12_data])

forward_solve = 0
chain_size = 1000
@time chain = sample(model, NUTS(0.65), chain_size,save_state=true)


#########################################
# The result from the sampled chain
#########################################
BurnIn = 0#00
chain_data = DataFrame(chain[BurnIn+1:end])
#################
# mean results
#################
freq_01_mean = mean(chain_data._freq_01)
freq_12_minus_mean = mean(chain_data._freq_12_minus)
freq_12_plus_mean = mean(chain_data._freq_12_plus)

fig_chain = plot(chain[BurnIn+1:end]);
xticks!(fig_chain[6],[3.44548800;3.44548880]);
xticks!(fig_chain[8],[3.23745050;3.23744950]);
display(fig_chain)
display(chain[BurnIn+1:end])


#########################################
# The result from the sampled chain
#########################################
chain_data = DataFrame(chain[BurnIn+1:end])
#################
# mean results
#################
freq_01_mean = mean(chain_data._freq_01)
freq_12_minus_mean = mean(chain_data._freq_12_minus)
freq_12_plus_mean  = mean(chain_data._freq_12_plus)

rho_chain_mean_01,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
                          2.0*pi*[freq_01_mean;freq_12_minus_mean],
                          2.0*pi*[freq_01_mean;freq_12_plus_mean],
                          [omr_R01,omegas[2]],
				   		  gamma1,gamma2,
						  0,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R01,
                          half_pi_option=half_pi_str)
population_chain_mean_01 = GLOQ.get_population(rho_chain_mean_01)
fig_mean_vs_noisy_01 = plot(t_R01./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_01,
		   		  label=["Mean-0" "Mean-1" "Mean-2"])
scatter!(fig_mean_vs_noisy_01,t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,
		 label=["Data-0" "Data-1" "Data-2"],legend=:outerright,
		 title = "Ramsey 0-1",
		 xlabel= "μs")
display(fig_mean_vs_noisy_01)

rho_chain_mean_12,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
                          2.0*pi*[freq_01_mean;freq_12_minus_mean],
                          2.0*pi*[freq_01_mean;freq_12_plus_mean],
                          [omegas[1];omr_R12],
				   		  gamma1,gamma2,
						  1,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)
population_chain_mean_12 = GLOQ.get_population(rho_chain_mean_12)
fig_mean_vs_noisy_12 = plot(t_R12./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_12,
		   		  label=["Mean-0" "Mean-1" "Mean-2"])
scatter!(fig_mean_vs_noisy_12,t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,
		 label=["Data-0" "Data-1" "Data-2"],legend=:outerright,
		 title = "Ramsey 1-2",
		 xlabel= "μs")
display(fig_mean_vs_noisy_12)


println( "Initial guess: freq 01 = ",freqs_opt[1],
		 " freq 12 minus = ",freqs_opt[2]-charge_noise_12,
         " freq 12 plus = ",freqs_opt[2]+charge_noise_12 )
println( "Mean values:   freq 01 = ",freq_01_mean,
		 " freq 12 minus = ",freqs_12_minus,
         " freq 12 plus = ",freqs_12_plus )


#########################
# Sample from the chain
#########################

global fig_result_01
global fig_result_12
num_random_sample = 200
for i = 1:num_random_sample
	global fig_result_01
	global fig_result_12
	sample_ind = rand(1:chain_size-BurnIn)
    freq_01_sample = chain_data._freq_01[sample_ind]
    freq_12_minus_sample = chain_data._freq_12_minus[sample_ind]
    freq_12_plus_sample = chain_data._freq_12_plus[sample_ind]

    rho_sample_u_01,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  (2.0*pi).*[freq_01_sample;freq_12_minus_sample],
                          (2.0*pi).*[freq_01_sample;freq_12_plus_sample],
                          [omr_R01;omegas[1]],
				   		  gamma1,gamma2,
						  0,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R01,
                          half_pi_option=half_pi_str)
    population_sample_01 = GLOQ.get_population(rho_sample_u_01)

	rho_sample_u_12,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  (2.0*pi).*[freq_01_sample;freq_12_minus_sample],
                          (2.0*pi).*[freq_01_sample;freq_12_plus_sample],
                          [omegas[1];omr_R12],
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
	  title=string("Ramsey 0-1, detuning = ",detuning_01)
	  )
plot!(fig_result_01,t_R01./GLOQ.GLOQ_MICRO_SEC,pop_R01_data,label="")
display(fig_result_01)

plot!(fig_result_12,t_R12./GLOQ.GLOQ_MICRO_SEC,population_chain_mean_12,label="",
	  xlabel="μs",
	  title=string("Ramsey 1-2, detuning = ",detuning_12)
		   		  )
plot!(fig_result_12,t_R12./GLOQ.GLOQ_MICRO_SEC,pop_R12_data,label="")
display(fig_result_12)

# Save the chain
date = string(today())

if (save_file)
	jls_file_name = string("chain-pm-",date,"/chain_",chain_size,".jls")
	write(jls_file_name,chain)
end

