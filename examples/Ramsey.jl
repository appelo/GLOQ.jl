using Random
using Plots
using DifferentialEquations,DiffEqFlux,LinearAlgebra
using Optim
using ForwardDiff
include("../src/GLOQ.jl")
pyplot()

N_states = 4;
omegas = (2.0*pi).*[4.10817777; 3.88303940; 3.61937188]
gamma1   = [2.59451982e-05; 4.54296537e-05; 0.0]
gamma2   = [2.58005604e-05; 9.00000000e-05; 0.0]
omr = 2.0*pi*(3.8830355 - 2.0e-3)
HK_free = GLOQ.RotationFrameDiagonal(omegas,omr)
L1,L2 = GLOQ.RotationFrameLindblad(gamma1,gamma2)

width = 17.00
TC = 2.5*width
Ω = 0.5*pi/(TC*sqrt(2.0))
θ = 0.0
HK_control,HS_control = GLOQ.RotationFrameRamseyControl(Ω,θ,N_states)



fromState = 1
rho_u0 = [0.0;0.0;0.0;0.0]
rho_v0 = [0.0;0.0;0.0;0.0]
rho_u0[fromState+1] = 1.0

T_Ramsey = 10.0 * 1000
N_dark_times = 101
t_dark_times = zeros(Float64,N_dark_times)
for i = 1:N_dark_times
	dark_time = (i-1)*dt
	t_dark_times[i] = dark_time
end
rho_Ramsey_u,rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr,
				 gamma1,gamma2,
				 TC,t_dark_times,N_states;initial_type = "states")
population_Ramsey = GLOQ.get_population(rho_Ramsey_u)


display(plot(t_dark_times./1000.0,population_Ramsey))
