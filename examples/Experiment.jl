using Random
using LinearAlgebra
using GalacticOptim,Optim,NLopt#,BlackBoxOptim
using DifferentialEquations#,DiffEqFlux
#using ForwardDiff
using Zygote
using Plots
#using ReverseDiff
include("../src/GLOQ.jl")
pyplot()
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

N_states = 4;
freqs = [4.10817777; 3.88303940; 3.61937188]
omegas = 2*pi.*freqs
gamma1   = [2.59451982e-05; 4.54296537e-05; 0.0]
gamma2   = [2.58005604e-05; 9.00000000e-05; 0.0]
omr_ramsey = 2.0*pi*(3.8830355 - 2.0e-3)
omr_echo = 2.0*pi*3.8830355
omr_T1 = 2.0*pi*3.8830355

width = 17.00
TC = 2.5*width
Ω = 0.5*pi/(TC*sqrt(2.0))
θ = 0.0

fromState = 1
rho_u0 = [0.0;0.0;0.0;0.0]
rho_v0 = [0.0;0.0;0.0;0.0]
rho_u0[fromState+1] = 1.0

T_Ramsey = 1.0*1000#10.0 * 1000
N_Ramsey = 11
t_ramsey = zeros(Float64,N_Ramsey)
dt_ramsey = T_Ramsey/(N_Ramsey-1)
for i = 1:N_Ramsey
	dark_time = (i-1)*dt_ramsey
	t_ramsey[i] = dark_time
end

# Ramsey Forward Solve
rho_ramsey_u_exp,rho_ramsey_v_exp = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr_ramsey,
				 gamma1,gamma2,
				 fromState,
				 TC,t_ramsey,N_states;
				 initial_type="states")
rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr_ramsey,
				 gamma1,gamma2,
				 fromState,
				 TC,t_ramsey,N_states;
				 method=STrapezoid())
population_Ramsey = GLOQ.get_population(rho_ramsey_u)
population_Ramsey_exp = GLOQ.get_population(rho_ramsey_u_exp)
fig=plot(t_ramsey./1000.0,population_Ramsey)
plot!(fig,t_ramsey./1000.0,population_Ramsey_exp)
display(fig)

fig2=plot(t_ramsey./1000.0,population_Ramsey-population_Ramsey_exp)
display(fig2)


#=
# Echo Forward Solve
T_Echo = 10.0*1000#10.0 * 1000
N_Echo = 1001
t_echo = zeros(Float64,N_Echo)
dt_echo = T_Echo/(N_Echo-1)
for i = 1:N_Echo
	dark_time = (i-1)*dt_echo
	t_echo[i] = dark_time
end
rho_echo_u_exp,rho_echo_v_exp = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
			     omegas,omr_echo,
				 gamma1,gamma2,
				 fromState,
				 TC,t_echo,N_states;
				 initial_type = "states")
rho_echo_u,rho_echo_v = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
			     omegas,omr_echo,
				 gamma1,gamma2,
				 fromState,
				 TC,t_echo,N_states;
				 initial_type = "states",
				 method = STrapezoid())
population_Echo = GLOQ.get_population(rho_echo_u)

fig2=plot(t_echo./1000.0,population_Echo)
display(fig2)

# T1 Forward Solve
T_T1 = 10.0*1000
N_T1 = 1001
t_T1 = zeros(Float64,N_T1)
dt_T1 = T_T1/(N_T1-1)
for i = 1:N_T1
	dark_time = (i-1)*dt_T1
	t_T1[i] = dark_time
end
rho_t1_u,rho_t1_v = GLOQ.T1ForwardSolve(rho_u0,rho_v0,
			     omegas,omr_T1,
				 gamma1,gamma2,
				 fromState,
				 TC,t_T1,N_states;
				 initial_type = "states")
population_T1 = GLOQ.get_population(rho_t1_u)

fig3=plot(t_T1./1000.0,population_T1)
display(fig3)
=#
