#using DiffEqSensitivity
using Zygote
#using ReverseDiff
using Random
using LinearAlgebra
#using DifferentialEquations#,DiffEqFlux
#using ForwardDiff
using Optim,NLopt
using Plots
#
using GLOQ
#include("../../src/GLOQ.jl")
pyplot()
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

N_states = 4;
freqs = [4.10817777; 3.88303940; 3.61937188]
omegas = 2.0*pi.*freqs
gamma1   = [2.59451982e-05; 4.54296537e-05; 0.0]
gamma2   = [2.58005604e-05; 9.00000000e-05; 0.0]
omr = 2.0*pi*(3.8830355)
HK_free = GLOQ.RotationFrameDiagonal(omegas,omr)
L1,L2 = GLOQ.RotationFrameLindblad(gamma1,gamma2)

width = 17.00
TC = 2.5*width

fromState = 1
rho_u0 = [0.0;0.0;0.0;0.0]
rho_v0 = [0.0;0.0;0.0;0.0]
rho_u0[fromState+1] = 1.0

T_Echo = 10.0*GLOQ.GLOQ_MICRO_SEC
N_dark_times = 26
t_dark_times = zeros(Float64,N_dark_times)
dt = T_Echo/(N_dark_times-1)
for i = 1:N_dark_times
	dark_time = (i-1)*dt
	t_dark_times[i] = dark_time
end
rho_Echo_u,rho_Echo_v = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
			     omegas,omr,
				 gamma1,gamma2,
				 1, # initial state
				 TC,t_dark_times,N_states)
population_Echo = GLOQ.get_population(rho_Echo_u)

fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_Echo)
display(fig)
#################################################
# Learn parameters
#################################################
p_true = [freqs;gamma1[1:2];gamma2[1:2]]
p_initial = [freqs.-1e-4;0.8.*gamma1[1:2];0.8.*gamma2[1:2]]
lower_bound = (0.7).*p_true
upper_bound = (1.3).*p_true

rho_Echo_u_guess,rho_Echo_v_guess = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
				 (2*pi).*p_initial[1:3],omr,
				 [p_initial[4:5];0.0],[p_initial[6:7];0.0],
				 1, # initial state
				 TC,t_dark_times,N_states;initial_type = "states")
population_guess= GLOQ.get_population(rho_Echo_u_guess)

plot!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_guess,line=(:dash),legend=:outerright,)
display(fig)
fig2 = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_guess-population_Echo,legend=:outerright,)
display(fig2)
##################################################################################################
global function_call
function loss(p)
	_rho_Echo_u,_rho_Echo_v = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
				     (2*pi).*p[1:3],omr,
					 [p[4:5];0.0],[p[6:7];0.0],#gamma1,gamma2,
					 1, # initial state
					 TC,t_dark_times,N_states;
					 )
	_population_Echo = GLOQ.get_population(_rho_Echo_u)
	_loss = sum(abs2,_population_Echo-population_Echo)/N_dark_times
	return _loss
end

function loss_gradient!(G,p)
	G .= Zygote.gradient(loss,p)[1]
end


println("Optim Fminbox(LBFGS) Optimization starts")
@time sol = Optim.optimize(loss,loss_gradient!,
								lower_bound,upper_bound,p_initial,
								Fminbox(LBFGS()),
								Optim.Options(
						        outer_iterations = 20,
								iterations = 10,
								show_trace=true,
								f_tol = 1e-3,
								outer_f_tol = 1e-3))
println("Optim Fminbox(LBFGS) Optimization done")


println("True p: ",p_true)
println("Initial p: ", p_initial," Initial error: ",p_initial-p_true)
println("\nOptim Fminbox-LBFGS: ",sol.minimizer," \nLoss:",sol.minimum," \nError: ",sol.minimizer-p_true)



#######################################################################

rho_Echo_u_optim,rho_Echo_v_optim = GLOQ.EchoForwardSolve(rho_u0,rho_v0,
				 (2*pi).*sol.minimizer[1:3],omr,
				 [sol.minimizer[4:5];0.0],[sol.minimizer[6:7];0.0],
				 1,# initial state
				 TC,t_dark_times,N_states;initial_type = "states")
population_optim = GLOQ.get_population(rho_Echo_u_optim)


plot!(fig2,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim-population_Echo,
	  line=(:dash),lab=:false,legend=:outerright,)
display(fig2)

fig3=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim-population_Echo,
	  line=(:dash),lab=:false,legend=:outerright,)
display(fig3)

fig4=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_Echo)
plot!(fig4,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim,line=(:dash))
