using Random
using Plots
using DifferentialEquations,DiffEqFlux,LinearAlgebra
using Optim,GalacticOptim
using ForwardDiff
using Zygote
include("../src/GLOQ.jl")
pyplot()

N_states = 4;
freqs = [4.10817777; 3.88303940; 3.61937188]
omegas = 2*pi.*freqs
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

T_Ramsey = 1.0*1000#10.0 * 1000
N_dark_times = 21#101
t_dark_times = zeros(Float64,N_dark_times)
dt = T_Ramsey/(N_dark_times-1)
for i = 1:N_dark_times
	dark_time = (i-1)*dt
	t_dark_times[i] = dark_time
end
rho_Ramsey_u,rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
			     omegas,omr,
				 gamma1,gamma2,
				 TC,t_dark_times,N_states;initial_type = "states")
population_Ramsey = GLOQ.get_population(rho_Ramsey_u)

fig=plot(t_dark_times./1000.0,population_Ramsey)
display(fig)

#################################################
# Learn parameters
#################################################
p_true = [freqs;gamma1[1:2];gamma2[1:2]]
p_initial = [freqs.-1e-4;0.9.*gamma1[1:2];0.9.*gamma2[1:2]]
lower_bound = (0.9).*p_true
upper_bound = (1.1).*p_true

rho_Ramsey_u_guess,rho_Ramsey_v_guess = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				 (2*pi).*p_initial[1:3],omr,
				 [p_initial[4:5];0.0],[p_initial[6:7];0.0],
				 TC,t_dark_times,N_states;initial_type = "states")
population_guess= GLOQ.get_population(rho_Ramsey_u_guess)

plot!(fig,t_dark_times./1000.0,population_guess,line=(:dash))
display(fig)
fig2 = plot(t_dark_times./1000.0,population_guess-population_Ramsey)
display(fig2)
##################################################################################################
function loss(p)
	_rho_Ramsey_u,_rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				     (2*pi).*p[1:3],omr,
					 [p[4:5];0.0],[p[6:7];0.0],#gamma1,gamma2,
					 TC,t_dark_times,N_states)
	_population_Ramsey = GLOQ.get_population(_rho_Ramsey_u)
	_loss = sum(abs2,_population_Ramsey-population_Ramsey)/N_dark_times
	return _loss
end
println("Initial loss: ",loss(p_initial))
@time p_DEFlux = DiffEqFlux.sciml_train(loss, p_initial, Fminbox(LBFGS()),
										lower_bounds=lower_bound,
										upper_bounds=upper_bound,
										show_trace=true,
										iterations=10,
										outer_iterations=10,
										f_tol=1e-4,
										outer_f_tol=5e-4,
										outer_x_tol=1e-6)


#@time Zygote.gradient(loss,p_initial)
#=
@time p_Optim = Optim.optimize(loss,
							   lower_bound,upper_bound,
							   p_initial,
							   Fminbox(LBFGS()),
							   Optim.Options(iterations=20,
											 x_tol = 1e-8,
											 g_tol = 1e-8,
											 outer_x_tol = 1e-6,
											 outer_g_tol = 1e-6,
							   				 show_trace=true)
							   )
=#

#=
adtype = GalacticOptim.AutoFiniteDiff()
@time p_DEFlux = DiffEqFlux.sciml_train(loss, p_initial, Fminbox(LBFGS()),
										adtype,
										lower_bounds=lower_bound,
										upper_bounds=upper_bound,
										show_trace=true,
										g_tol=1e-6,
										x_tol=1e-6,
										iterations=25)
=#

println("True p: ",p_true)
println("Initial p: ", p_initial," Initial error: ",p_initial-p_true)
println("DEFlux p: ",p_DEFlux.u," DEFlux error: ",p_DEFlux.u-p_true)
println(p_DEFlux)


#######################################################################
rho_Ramsey_u_DEFlux,rho_Ramsey_v_DEFlux = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				 (2*pi).*p_DEFlux.u[1:3],omr,
				 [p_DEFlux.u[4:5];0.0],[p_DEFlux.u[6:7];0.0],
				 TC,t_dark_times,N_states;initial_type = "states")
population_DEFlux = GLOQ.get_population(rho_Ramsey_u_DEFlux)


plot!(fig2,t_dark_times./1000.0,population_DEFlux-population_Ramsey,line=(:dash),lab=:false)
display(fig2)

fig3=plot(t_dark_times./1000.0,population_Ramsey)
plot!(fig3,t_dark_times./1000.0,population_DEFlux,line=(:dash))


#println("Optim p: ",p_Optim.minimizer)
#println(p_Optim)
