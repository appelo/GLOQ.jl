using DiffEqSensitivity
using Zygote
using Random
using LinearAlgebra
using DifferentialEquations
using NLopt,Optim
using Plots
#include("../../src/GLOQ.jl")
using GLOQ
pyplot()
# For this example, NLopt's LBFGS with box constrains seems to be faster than Optim's FMinbox(LBFGS)

N_states = 4;
freqs = [4.10817777; 3.88303940; 3.61937188]
omegas = 2.0*pi.*freqs
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

T_Ramsey = 1.0*GLOQ.GLOQ_MICRO_SEC
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
				 1, # initial state
				 TC,t_dark_times,N_states)
population_Ramsey = GLOQ.get_population(rho_Ramsey_u)

fig=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_Ramsey)
#display(fig)
#################################################
# Learn parameters
#################################################
p_true = [freqs;gamma1[1:2];gamma2[1:2]]
p_initial = [freqs.-1e-4;0.9.*gamma1[1:2];0.9.*gamma2[1:2]]
lower_bound = (0.85).*p_true
upper_bound = (1.15).*p_true

rho_Ramsey_u_guess,rho_Ramsey_v_guess = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				 (2*pi).*p_initial[1:3],omr,
				 [p_initial[4:5];0.0],[p_initial[6:7];0.0],
				 1, # initial state
				 TC,t_dark_times,N_states;initial_type = "states") 
population_guess= GLOQ.get_population(rho_Ramsey_u_guess)

plot!(fig,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_guess,line=(:dash),legend=:outerright,)
display(fig)
fig2 = plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_guess-population_Ramsey,legend=:outerright,)
display(fig2)
##################################################################################################
function loss_optim(p)
	_rho_Ramsey_u,_rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				     (2*pi).*p[1:3],omr,
					 [p[4:5];0.0],[p[6:7];0.0],#gamma1,gamma2,
					 1, # initial state
					 TC,t_dark_times,N_states;#)
					 initial_type="states")
	_population_Ramsey = GLOQ.get_population(_rho_Ramsey_u)
	_loss = sum(abs2,_population_Ramsey-population_Ramsey)/N_dark_times
	return _loss
end

function gradient_optim!(G,x)
	G .= Zygote.gradient(loss_optim,x)[1]
end

function loss_nlopt(x,grad)
	grad .= Zygote.gradient(loss_optim,x)[1]
	return loss_optim(x)
end

println("Optim Fminbox(LBFGS) Optimization starts")
@time sol_optim_minbox_lbfgs = Optim.optimize(loss_optim,gradient_optim!,
								lower_bound,upper_bound,p_initial,
								Fminbox(LBFGS()),
						        Optim.Options(outer_iterations = 20,
								iterations = 10,
								show_trace=false,
								f_tol = 1e-3,outer_f_tol = 1e-3))
println("Optim Fminbox(LBFGS) Optimization done")


function_call_num = 0
println("NLopt LBFGS Optimization with Zygote gradient starts")
nl_opt_obj = NLopt.Opt(:LD_LBFGS,length(p_initial))
nl_opt_obj.min_objective = loss_nlopt
nl_opt_obj.lower_bounds = lower_bound
nl_opt_obj.upper_bounds = upper_bound
nl_opt_obj.ftol_rel = 1e-3
#nl_opt_obj.maxeval = 200
@time sol_nlopt_LBFGS = NLopt.optimize(nl_opt_obj,p_initial)
println("NLopt LBFGS Optimization with Zygote gradient done")


println("True p: ",p_true)
println("Initial p: ", p_initial," Initial error: ",p_initial-p_true)
#println("\nOptim Fminbox-LBFGS: ",sol_optim_minbox_lbfgs.minimizer," \nLoss:",sol_optim_minbox_lbfgs.minimum," \nError: ",sol_optim_minbox_lbfgs.minimizer-p_true)
println("\nNLopt LBFGS: ",sol_nlopt_LBFGS[2]," \nLoss:",sol_nlopt_LBFGS[1]," \nError: ",sol_nlopt_LBFGS[2]-p_true)
#println("\nNLopt LBFGS with ReverseDiff gradient: ",sol_nlopt_LBFGS_rd.u," \nLoss:",sol_nlopt_LBFGS_rd.minimum," \nError: ",sol_nlopt_LBFGS_rd.u-p_true)




#######################################################################

rho_Ramsey_u_optim,rho_Ramsey_v_optim = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				 (2*pi).*sol_optim_minbox_lbfgs.minimizer[1:3],omr,
				 [sol_optim_minbox_lbfgs.minimizer[4:5];0.0],[sol_optim_minbox_lbfgs.minimizer[6:7];0.0],
				 1,# initial state
				 TC,t_dark_times,N_states;initial_type = "states")
population_optim = GLOQ.get_population(rho_Ramsey_u_optim)


plot!(fig2,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim-population_Ramsey,
	  line=(:dash),lab=:false,legend=:outerright)
display(fig2)

fig3=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim-population_Ramsey,
	  line=(:dash),lab=:false,legend=:outerright)
display(fig3)

fig4=plot(t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_Ramsey)
plot!(fig4,t_dark_times./GLOQ.GLOQ_MICRO_SEC,population_optim,line=(:dash))
