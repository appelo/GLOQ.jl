using Random
using LinearAlgebra
using Optim,NLopt
using DifferentialEquations
using Zygote
using Plots
using GLOQ
pyplot()

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

T_Ramsey = 1.0*GLOQ.GLOQ_MICRO_SEC#10.0 * 1000
N_Ramsey = 101
t_ramsey = zeros(Float64,N_Ramsey)
dt_ramsey = T_Ramsey/(N_Ramsey-1)
for i = 1:N_Ramsey
	dark_time = (i-1)*dt_ramsey
	t_ramsey[i] = dark_time
end

# Ramsey Forward Solve
@time rho_ramsey_u_exp,rho_ramsey_v_exp = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr_ramsey,
				 gamma1,gamma2,
				 fromState,
				 TC,t_ramsey,N_states;
				 initial_type="states")
@time rho_ramsey_u,rho_ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr_ramsey,
				 gamma1,gamma2,
				 fromState,
				 TC,t_ramsey,N_states;
				 method = Trapezoid(),
				 abstol = 1e-10,
				 )
population_Ramsey = GLOQ.get_population(rho_ramsey_u)
population_Ramsey_exp = GLOQ.get_population(rho_ramsey_u_exp)

println("Ramsey max error:",maximum(abs.(population_Ramsey-population_Ramsey_exp)) )

fig=plot(t_ramsey./GLOQ.GLOQ_MICRO_SEC,population_Ramsey)
plot!(fig,t_ramsey./GLOQ.GLOQ_MICRO_SEC,population_Ramsey_exp)
display(fig)

fig2=plot(t_ramsey./GLOQ.GLOQ_MICRO_SEC,population_Ramsey-population_Ramsey_exp)
display(fig2)



# Echo Forward Solve
T_Echo = 10.0*GLOQ.GLOQ_MICRO_SEC#10.0 * 1000
N_Echo = 1001
t_echo = zeros(Float64,N_Echo)
dt_echo = T_Echo/(N_Echo-1)
for i = 1:N_Echo
	dark_time = (i-1)*dt_echo
	t_echo[i] = dark_time
end
@time rho_echo_u_exp,rho_echo_v_exp = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
			     omegas,omr_echo,
				 gamma1,gamma2,
				 fromState,
				 TC,t_echo,N_states;
				 initial_type = "states")
@time rho_echo_u,rho_echo_v = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
			     omegas,omr_echo,
				 gamma1,gamma2,
				 fromState,
				 TC,t_echo,N_states;
				 initial_type = "states",
				 method = Trapezoid(),
				 abstol = 1e-10)
population_Echo = GLOQ.get_population(rho_echo_u)
population_Echo_exp = GLOQ.get_population(rho_echo_u_exp)

println( maximum(abs.(population_Echo-population_Echo_exp)) )

fig_echo=plot(t_echo./GLOQ.GLOQ_MICRO_SEC,population_Echo)
plot!(fig_echo,t_echo./GLOQ.GLOQ_MICRO_SEC,population_Echo_exp,line=(:dash))
display(fig_echo)

fig_echo2=plot(t_echo./GLOQ.GLOQ_MICRO_SEC,population_Echo-population_Echo_exp)
display(fig_echo2)

# T1 Forward Solve
T_T1 = 10.0*1000
N_T1 = 1001
t_T1 = zeros(Float64,N_T1)
dt_T1 = T_T1/(N_T1-1)
for i = 1:N_T1
	dark_time = (i-1)*dt_T1
	t_T1[i] = dark_time
end
@time rho_t1_u_exp,rho_t1_v_exp = GLOQ.T1ForwardSolve(rho_u0,rho_v0,
			     omegas,omr_T1,
				 gamma1,gamma2,
				 fromState,
				 TC,t_T1,N_states;
				 initial_type = "states")
population_T1_exp = GLOQ.get_population(rho_t1_u_exp)
@time rho_t1_u,rho_t1_v = GLOQ.T1ForwardSolve(rho_u0,rho_v0,
			     omegas,omr_T1,
				 gamma1,gamma2,
				 fromState,
				 TC,t_T1,N_states;
				 initial_type = "states",
				 method = Trapezoid(),
				 abstol = 1e-10)
population_T1 = GLOQ.get_population(rho_t1_u)

println( maximum(abs.(population_T1-population_T1_exp)) )

fig_t1=plot(t_T1./GLOQ.GLOQ_MICRO_SEC,population_T1)
plot!(fig_t1,t_T1./GLOQ.GLOQ_MICRO_SEC,population_T1_exp,line=(:dash))
display(fig_t1)

fig_t1_error=plot(t_T1./GLOQ.GLOQ_MICRO_SEC,population_T1-population_T1_exp)
display(fig_t1_error)
