

#############################################################################
# Read data files
#############################################################################
half_pi_str = "duration" 
data_set_option = "short"
#data_set_option = "long"

# Load Ramsey 01 data
if (data_set_option=="short")
#    ramsey_01_0_file = "/Users/zhichaopeng/Dropbox/research_projects/Quantum/Experiment/UQ/0413/0413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_0.txt"
    ramsey_01_0_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_0.txt"
    ramsey_01_1_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_1.txt"
    ramsey_01_2_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_4000_5_2.txt"
    t_ramsey_01_file = "data-set-20220413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_4000_5.txt"

    ramsey_12_0_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_0.txt"
    ramsey_12_1_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_1.txt"
    ramsey_12_2_file = "data-set-20220413/population_ramsey_12_1000000.0_4000_5_1000_confusion_2.txt"
    t_ramsey_12_file = "data-set-20220413/darktime_ramsey_12_1000000.0_5_4000_1000_20220413_confusion.txt"
else
    ramsey_01_0_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_0.txt"
    ramsey_01_1_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_1.txt"
    ramsey_01_2_file = "data-set-20220413/population_ramsey_01_half_period_1000000.0_1000_20220413_80000_20_2.txt"
    t_ramsey_01_file = "data-set-20220413/darktime_ramsey_01_half_period_1000000.0_1000_20220413_80000_20.txt"

    ramsey_12_0_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_0.txt"
    ramsey_12_1_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_1.txt"
    ramsey_12_2_file = "data-set-20220413/population_ramsey_12_1000000.0_80000_20_1000_confusion_2.txt"
    t_ramsey_12_file = "data-set-20220413/darktime_ramsey_12_1000000.0_20_80000_1000_20220413_confusion.txt"
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
omegas = 2.0*pi.*freqs # change to angular frequency

freqs_opt = [3.445488401954036;3.2374501243705796] # transition frequency charaterized with deterministic code
omegas_opt = 2.0*pi.*freqs_opt

T1 = [169.1;85.79]
gamma1   = [1.0/(T1[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T1[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of relaxation time - T1 (in units of ns)
#T2 = [7.5;2.5]
T2 = [6.9309771664384865;1.965597856818199]
charge_noise_01 = 0.0
charge_noise_12 = 1.229132624029905e-4
gamma2   = [1.0/(T2[1]*GLOQ.GLOQ_MICRO_SEC); 1.0/(T2[2]*GLOQ.GLOQ_MICRO_SEC)] # Reciprocal of dephasing time - T2 (in units of ns)

charge_noise = 2.0*pi*[charge_noise_01;charge_noise_12]

freqs_12_minus = freqs_opt[2]-charge_noise_12
freqs_12_plus  = freqs_opt[2]+charge_noise_12
omegas_minus = omegas_opt-charge_noise
omegas_plus  = omegas_opt+charge_noise
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

#######################################################
# Set up the autodifferentiation back end for Turing
#######################################################
Turing.setadbackend(:zygote)
#Turing.setadbackend(:reversediff)
#Turing.setadbackend(:forwarddiff)
global forward_solve
forward_solve = 0

N_dark_R01 = length(pop_R01_0)
N_dark_R12 = length(pop_R12_0)

@model function RamseyExperiment_01_12(data)
	 # Priori distribution
	#σ ~ Normal(0.0,0.5)
	σ_01 ~ InverseGamma()
	σ_12 ~ InverseGamma()
    # obtained from the determinisitic optimization
    _freq_01 ~ truncated(Normal(freqs_opt[1],5e-5),freqs_opt[1]-5e-4,freqs_opt[1]+5e-4)
	_freq_12_minus ~ truncated(Normal(freqs_12_minus,5e-5),freqs_12_minus-1.25e-4,freqs_12_minus+1.25e-4)
    _freq_12_plus  ~ truncated(Normal(freqs_12_plus,5e-5), freqs_12_plus-1.25e-4,freqs_12_plus+1.25e-4)

    _omegas_minus = (2.0*pi).*[_freq_01;_freq_12_minus]
    _omegas_plus  = (2.0*pi).*[_freq_01;_freq_12_plus]

	_rho_ramsey_01,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  _omegas_minus,_omegas_plus,
                          [omr_R01;omegas[2]],
				   		  gamma1,gamma2,
						  0,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R01,
                          half_pi_option=half_pi_str)

	_rho_ramsey_12,_ = GLOQ.PerformForwardRamseyExperimentParityPlusMinus(
                          state_u0,state_v0,
				   		  _omegas_minus,_omegas_plus,
                          [omegas[1];omr_R12],
				   		  gamma1,gamma2,
						  1,
				   		  [TC01;TC12],
						  [amp_01;amp_12],
						  t_R12,
                          half_pi_option=half_pi_str)

	_population_ramsey_01 = GLOQ.get_population(_rho_ramsey_01)
	_population_ramsey_12 = GLOQ.get_population(_rho_ramsey_12)

	for i = 1:N_dark_R01 
        data[i,:] ~ MvNormal(_population_ramsey_01[i,:], σ_01)
    end

	for i = 1:N_dark_R12
		data[N_dark_R01+i,:] ~ MvNormal(_population_ramsey_12[i,:], σ_12)
	end
	global forward_solve
	forward_solve += 1
	println("Forward solve ",forward_solve," done")
end
