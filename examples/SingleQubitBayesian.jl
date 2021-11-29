using LinearAlgebra
using Zygote,ReverseDiff
using Turing, Distributions, DifferentialEquations
# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots, StatsPlots
using DataFrames
# Set a seed for reproducibility.
using Random
Random.seed!(14);
using Plots
include("../src/GLOQ.jl")
pyplot()

# System parameters for a simple two level open quantum system
N_states = 2; # number of states
freqs = [4.1] # transition frequency in GHz
omegas = 2.0*pi.*freqs # change to angular frequency
gamma1   = [25e-05] # decay???
gamma2   = [25e-05] # dephasing???
omr = 2.0*pi*(4.1 - 1.0e-3) # drive frequency
TC = 2.5*17.0 # total control time

# Initial state
initial_state = 0
rho_u0 = [0.0;0.0]
rho_v0 = [0.0;0.0]
rho_u0[initial_state+1] = 1.0

# Duration of the Ramsey experiment, largest dark time
T_Ramsey = 10.0*GLOQ.GLOQ_MICRO_SEC # convert micro-sec to nano-sec
# total number of dark time samples
N_dark_times = 201
t_dark_times = collect(range(0.0, T_Ramsey, length=N_dark_times))

# Forward solve to generate synthetic data
rho_synthetic_ramsey_u,rho_synthetic_ramsey_v = GLOQ.RamseyForwardSolve(
				 rho_u0,rho_v0, # initial values, u for the real part, v for the imaginary part
			     omegas,omr, # transition frequencies, drive frequency
				 gamma1,gamma2, # decay and dephasing parameters ?
				 initial_state, # initial state
				 TC,t_dark_times,N_states) # control time, dark time, total number of states
population_synthetic = GLOQ.get_population(rho_synthetic_ramsey_u)

# Add noise to the synthetic data
multiplicative_noise = 1.0.+0.025*randn(N_dark_times)
additive_noise = 0.05*randn(N_dark_times)
noisy_data = population_synthetic.*multiplicative_noise
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

println(size(noisy_data))
# plot the synthetic data with and without noise
fig = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic)
scatter!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,noisy_data)
display(fig)

#fig2 = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic-noisy_data)
#display(fig2)

#noisy_data = transpose(noisy_data)
p_true = [freqs;gamma1;gamma2]

#######################################################
# Set up the autodifferentiation back end for Turing
#Turing.setadbackend(:zygote)
Turing.setadbackend(:reversediff)
#Turing.setadbackend(:forwarddiff)
global sample_number
@model function RamseyExperiment(data)
	 σ ~ InverseGamma()
    _freq ~ truncated(Normal(4.1,2e-4),4.1-5e-4,4.1+5e-4)
    _gamma1 ~ truncated(Normal(25e-05,1e-5),10e-5,40e-5)
    _gamma2 ~ truncated(Normal(25e-05,1e-5),10e-5,40e-5)

	_rho_ramsey_u,_rho_ramsey_v = GLOQ.RamseyForwardSolve(
					 rho_u0,rho_v0, # initial values, u for the real part, v for the imaginary part
				     2.0*pi*[_freq],omr, # transition frequencies, drive frequency
					 [_gamma1],[_gamma2], # decay and dephasing parameters ?
					 initial_state, # initial state
					 TC,t_dark_times,N_states;
					 method="exponential")
					 #method = Trapezoid()) # control time, dark time, total number of states
	_population_ramsey = GLOQ.get_population(rho_synthetic_ramsey_u)
	for i = 1:N_dark_times
        data[i,:] ~ MvNormal(_population_ramsey[i,:], σ)
    end
	global sample_number
	sample_number += 1
	println("Sample ",sample_number," done")
end

model = RamseyExperiment(noisy_data)

sample_number = 0
@time chain = sample(model, MH(),500)
display(plot(chain))
display(chain)

# extact the data value from the chain
chain_data = DataFrame(chain)
freq_mean = mean(chain_data._freq)
gamma1_mean = mean(chain_data._gamma1)
gamma2_mean = mean(chain_data._gamma2)
rho_chain_mean_u,rho_chain_mean_v = GLOQ.RamseyForwardSolve(
				 rho_u0,rho_v0, # initial values, u for the real part, v for the imaginary part
				 2.0*pi*[freq_mean],omr, # transition frequencies, drive frequency
				 [gamma1_mean],[gamma2_mean], # decay and dephasing parameters ?
				 initial_state, # initial state
				 TC,t_dark_times,N_states;
				 method="exponential")
population_chain_mean = GLOQ.get_population(rho_chain_mean_u)
fig_result = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_chain_mean,label=["Mean-0" "Mean-1"])
scatter!(fig_result,t_dark_times./GLOQ.GLOQ_MICRO_SEC,noisy_data,label=["Noisy data-0" "Noisy data-1"])
display(fig_result)

fig_result2 = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_chain_mean,label=["Mean-0" "Mean-1"])
scatter!(fig_result2,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_synthetic,label=["Data-0" "Data-1"])
display(fig_result2)
