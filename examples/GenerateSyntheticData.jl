using Zygote
using Random
using DelimitedFiles
using LinearAlgebra
using GalacticOptim,Optim,NLopt#,BlackBoxOptim
using DifferentialEquations#,DiffEqFlux
#using ForwardDiff
using Plots
using LaTeXStrings
#using GLOQ
#using ReverseDiff
include("../src/GLOQ.jl")
pyplot()
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

N_states = 4;
freqs = [4.0108; 3.8830; 3.6287]
omegas = 2.0*pi.*freqs
gamma1   = [2.222222e-05; 4.761904e-05; 4.545455e-5]
gamma2   = [4.166667e-05; 6.8965517e-5; 2.3255814e-4]

# driving frequency
omr_t1_01 = 2.0*pi*4.0108
omr_t1_12 = 2.0*pi*3.8830
omr_t1_23 = 2.0*pi*3.6287

omr_ramsey_01 = 2.0*pi*(4.0108-5e-4)
omr_ramsey_12 = 2.0*pi*(3.8830-2.0e-3)
omr_ramsey_23 = 2.0*pi*(3.6287-2.0e-3)

omr_echo_01 = 2.0*pi*4.0108
omr_echo_12 = 2.0*pi*3.8830
omr_echo_23 = 2.0*pi*3.6287

# charge noise
char_noise = 2.0*pi.*[2e-6,6e-5,1.3e-4]
# width of the control signal
width = 21.25
TC_01 = 2.5*width

width = 17.00
TC_12 = 2.5*width

width = 17.00
TC_23 = 2.5*width

# in nano second
T_T1_01 = 125.0*1000.0
T_T1_12 = 125.0*1000.0
T_T1_23 = 125.0*1000.0
T_Ramsey_01 = 20.0*1000.0
T_Ramsey_12 = 20.0*1000.0
T_Ramsey_23 = 10.0*1000.0
T_Echo_01 = 40.0*1000.0
T_Echo_12 = 20.0*1000.0
T_Echo_23 = 10.0*1000.0

# data size used for characterization
n1use = 201
n2use = 201
n3use = 201
n4yse = 201
n4use = 201
n5use = 201
n6use = 201
n7use = 201
n8use = 201
n9use = 201

# time points evaluated
dt_t1_01 = T_T1_01/200
dt_t1_12 = T_T1_12/200
dt_t1_23 = T_T1_23/200
dt_ramsey_01 = T_Ramsey_01/200
dt_ramsey_12 = T_Ramsey_12/200
dt_ramsey_23= T_Ramsey_23/200
dt_echo_01 = T_Echo_01/200
dt_echo_12 = T_Echo_12/200
dt_echo_23 = T_Echo_23/200

t_t1_01 = zeros(Float64,n1use)
for i = 1:n1use
	t_t1_01[i] = (i-1)*dt_t1_01
end
T_evaluated_t1_01 = t_t1_01[end]

t_t1_12 = zeros(Float64,n2use)
for i = 1:n2use
	t_t1_12[i] = (i-1)*dt_t1_12
end
T_evaluated_t1_12 = t_t1_12[end]

t_t1_23 = zeros(Float64,n3use)
for i = 1:n3use
	t_t1_23[i] = (i-1)*dt_t1_23
end
T_evaluated_t1_23 = t_t1_23[end]

t_ramsey_01 = zeros(Float64,n4use)
for i = 1:n4use
	t_ramsey_01[i] = (i-1)*dt_ramsey_01
end
T_evaluated_ramsey_01 = t_ramsey_01[end]

t_ramsey_12 = zeros(Float64,n5use)
for i = 1:n5use
	t_ramsey_12[i] = (i-1)*dt_ramsey_12
end
T_evaluated_ramsey_12 = t_ramsey_12[end]

t_ramsey_23 = zeros(Float64,n6use)
for i = 1:n6use
	t_ramsey_23[i] = (i-1)*dt_ramsey_23
end
T_evaluated_ramsey_23 = t_ramsey_23[end]

t_echo_01 = zeros(Float64,n7use)
for i = 1:n7use
	t_echo_01[i] = (i-1)*dt_echo_01
end
T_evaluated_echo_01 = t_echo_01[end]

t_echo_12 = zeros(Float64,n8use)
for i = 1:n8use
	t_echo_12[i] = (i-1)*dt_echo_12
end
T_evaluated_echo_12 = t_echo_12[end]

t_echo_23 = zeros(Float64,n9use)
for i = 1:n9use
	t_echo_23[i] = (i-1)*dt_echo_23
end
T_evaluated_echo_23 = t_echo_23[end]


# initial conditions
u0_t1_01 = [1.0;0.0;0.0;0.0]
u0_t1_12 = [0.0;1.0;0.0;0.0]
u0_t1_23 = [0.0;0.0;1.0;0.0]
u0_ramsey_01 = [1.0;0.0;0.0;0.0]
u0_ramsey_12 = [0.0;1.0;0.0;0.0]
u0_ramsey_23 = [0.0;0.0;1.0;0.0]
u0_echo_01 = [1.0;0.0;0.0;0.0]
u0_echo_12 = [0.0;1.0;0.0;0.0]
u0_echo_23 = [0.0;0.0;1.0;0.0]

v0_t1_01 = [0.0;0.0;0.0;0.0]
v0_t1_12 = [0.0;0.0;0.0;0.0]
v0_t1_23 = [0.0;0.0;0.0;0.0]
v0_ramsey_01 = [0.0;0.0;0.0;0.0]
v0_ramsey_12 = [0.0;0.0;0.0;0.0]
v0_ramsey_23 = [0.0;0.0;0.0;0.0]
v0_echo_01 = [0.0;0.0;0.0;0.0]
v0_echo_12 = [0.0;0.0;0.0;0.0]
v0_echo_23 = [0.0;0.0;0.0;0.0]

########################################################
# ForwardSolves
########################################################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 omegas,omr_echo_01,
				 char_noise,
				 gamma1,gamma2,
				 0,
				 TC_01,t_echo_01,N_states)
population_echo_01 = GLOQ.get_population(rho_Echo01_u)

# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 omegas,omr_echo_12,
				 char_noise,
				 gamma1,gamma2,
				 1,
				 TC_12,t_echo_12,N_states)
population_echo_12 = GLOQ.get_population(rho_Echo12_u)

# Echo 2-3
rho_Echo23_u,rho_Echo23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
				 omegas,omr_echo_23,
				 char_noise,
				 gamma1,gamma2,
				 2,
				 TC_23,t_echo_23,N_states)
population_echo_23 = GLOQ.get_population(rho_Echo23_u)

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 omegas,omr_ramsey_01,
				 char_noise,
				 gamma1,gamma2,
				 0,
				 TC_01,t_ramsey_01,N_states)
population_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)

# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 omegas,omr_ramsey_12,
				 char_noise,
				 gamma1,gamma2,
				 1,
				 TC_12,t_ramsey_12,N_states)
population_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)

# Ramsey 2-3
@time rho_Ramsey23_u,rho_Ramsey23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
				 omegas,omr_ramsey_23,
				 char_noise,
				 gamma1,gamma2,
				 2,
				 TC_23,t_ramsey_23,N_states)
population_ramsey_23 = GLOQ.get_population(rho_Ramsey23_u)

# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 omegas,omr_t1_01,
				 char_noise,
				 gamma1,gamma2,
				 0,
				 TC_01,t_t1_01,N_states)
population_t1_01 = GLOQ.get_population(rho_T101_u)

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 omegas,omr_t1_12,
				 char_noise,
				 gamma1,gamma2,
				 1,
				 TC_12,t_t1_12,N_states)
population_t1_12 = GLOQ.get_population(rho_T112_u)

# T1 2-3
rho_T123_u,rho_T123_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
				 omegas,omr_t1_23,
				 char_noise,
				 gamma1,gamma2,
				 2,
				 TC_23,t_t1_23,N_states)
population_t1_23 = GLOQ.get_population(rho_T123_u)

# Plotting
fnt = Plots.font("Helvetica",16)
lfnt = Plots.font("Helvetica",12)
Plots.default(titlefont=fnt,
			  guidefont=fnt,
			  tickfont=fnt,
			  legendfont=lfnt,
			  linewidth=2);

# Echo 0-1
fig_e01 = plot(t_echo_01./1000.0,population_echo_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_e01)


# Echo 1-2
fig_e12 = plot(t_echo_12./1000.0,population_echo_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_e12)

# Echo 2-3
fig_e23 = plot(t_echo_23./1000.0,population_echo_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_e23)


# Ramsey 0-1
fig_r01 = plot(t_ramsey_01./1000.0,population_ramsey_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_r01)


# Ramsey 1-2
fig_r12 = plot(t_ramsey_12./1000.0,population_ramsey_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_r12)

# Ramsey 2-3
fig_r23 = plot(t_ramsey_23./1000.0,population_ramsey_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_r23)

# T1 0-1
fig_t01 = plot(t_t1_01./1000.0,population_t1_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_t01)

# T1 1-2
fig_t12 = plot(t_t1_12./1000.0,population_t1_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_t12)

# T1 1-2
fig_t23 = plot(t_t1_23./1000.0,population_t1_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
display(fig_t23)

# Save data
writedlm( string("synthetic_data/echo_01.txt"),  population_echo_01, ',')
writedlm( string("synthetic_data/echo_12.txt"),  population_echo_12, ',')
writedlm( string("synthetic_data/echo_23.txt"),  population_echo_23, ',')

#writedlm( string("synthetic_data/ramsey_01.txt"),  population_ramsey_01, ',')
#writedlm( string("synthetic_data/ramsey_12.txt"),  population_ramsey_12, ',')
#writedlm( string("synthetic_data/ramsey_23.txt"),  population_ramsey_23, ',')

writedlm( string("synthetic_data/t1_01.txt"),  population_t1_01, ',')
writedlm( string("synthetic_data/t1_12.txt"),  population_t1_12, ',')
writedlm( string("synthetic_data/t1_23.txt"),  population_t1_23, ',')
