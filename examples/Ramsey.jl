#using Zygote
using ReverseDiff
using Random
using LinearAlgebra
using DifferentialEquations#,DiffEqFlux
using DiffEqSensitivity
#using ForwardDiff
using GalacticOptim,Optim,NLopt
using Plots
#
#using GLOQ
include("../src/GLOQ.jl")
pyplot()
# BlackBoxOptim somehow downgrade some packages and as a result breaks the auto-differentiation with Zygote
# we should avoid it.
# NLopt seems to be really fast. Their LBFGS with box constrains is faster than Optim's FMinbox(LBFGS)

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
				 1, # initial state
				 TC,t_dark_times,N_states)
population_Ramsey = GLOQ.get_population(rho_Ramsey_u)

fig=plot(t_dark_times./1000.0,population_Ramsey)
display(fig)
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

plot!(fig,t_dark_times./1000.0,population_guess,line=(:dash),legend=:outerright,)
display(fig)
fig2 = plot(t_dark_times./1000.0,population_guess-population_Ramsey,legend=:outerright,)
display(fig2)
##################################################################################################
global function_call
function loss(p)
	_rho_Ramsey_u,_rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				     (2*pi).*p[1:3],omr,
					 [p[4:5];0.0],[p[6:7];0.0],#gamma1,gamma2,
					 1, # initial state
					 TC,t_dark_times,N_states;
					 #method=Tsit5() )
					 method="DiffEqDefault")
	_population_Ramsey = GLOQ.get_population(_rho_Ramsey_u)
	#println(size(_population_Ramsey))
	#println(size(population_Ramsey))
	_loss = sum(abs2,_population_Ramsey-population_Ramsey)/N_dark_times
	return _loss
end
println("Initial loss: ",loss(p_initial))
@time ReverseDiff.gradient(loss,p_initial)

global function_call_num
function_call_num = 0
function loss_gala(p,pp)
	global function_call_num
	function_call_num += 1
	_rho_Ramsey_u,_rho_Ramsey_v = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				     (2*pi).*p[1:3],omr,
					 [p[4:5];0.0],[p[6:7];0.0],#gamma1,gamma2,
					 1, # initial state
					 TC,t_dark_times,N_states;
					 #method = "exponential")
					 #method="DiffEqDefault")
					 method=STrapezoid())
	_population_Ramsey = GLOQ.get_population(_rho_Ramsey_u)
	_loss = sum(abs2,_population_Ramsey-population_Ramsey)/N_dark_times
	#println("Function call ",function_call_num," Loss: ",_loss)
	return _loss
end

#=
struct OptimizationSolution{T, N, uType, P, A, Tf, O} <: AbstractOptimizationSolution{T, N}
    u::uType # minimizer
    prob::P # optimization problem
    alg::A # algorithm
    minimum::Tf
    retcode::Symbol
    original::O # original output of the optimizer
end
=#

optimization_function_fd = OptimizationFunction(loss_gala, GalacticOptim.AutoFiniteDiff())
prob_fd = GalacticOptim.OptimizationProblem(optimization_function_fd, p_initial,
										 lb = lower_bound, ub = upper_bound)

optimization_function_autozygote = OptimizationFunction(loss_gala, GalacticOptim.AutoZygote())
prob_autozygote = GalacticOptim.OptimizationProblem(optimization_function_autozygote, p_initial,
										 lb = lower_bound, ub = upper_bound)

function_call_num = 0

println("Optim Fminbox(LBFGS) Optimization starts")
@time sol_optim_minbox_lbfgs = GalacticOptim.solve(prob_autozygote ,Fminbox(LBFGS()),
						        outer_iterations = 20,
								iterations = 10,
								show_trace=true,
								f_tol = 1e-3,
								outer_f_tol = 1e-3)
println("Optim Fminbox(LBFGS) Optimization done")

function_call_num = 0
println("NLopt Nelder-Mead Optimization starts")
@time sol_nlopt_NedlerMead = solve(prob_autozygote , Opt(:LN_NELDERMEAD,length(p_initial)),
								   maxiters=200,
								   ftol_rel=1e-3)
println("NLopt Nelder-Mead Optimization done")

function_call_num = 0
println("NLopt LBFGS Optimization with Zygote gradient starts")
@time sol_nlopt_LBFGS = solve(prob_autozygote , Opt(:LD_LBFGS,length(p_initial)),
								   maxiters=200,
								   ftol_rel=1e-3)
println("NLopt LBFGS Optimization with Zygote gradient done")

function_call_num = 0
println("NLopt LBFGS Optimization with FD gradient starts")
@time sol_nlopt_LBFGS_fd = solve(prob_fd , Opt(:LD_LBFGS,length(p_initial)),
								   maxiters=200,
								   ftol_rel=1e-3)
println("NLopt LBFGS Optimization with FD gradient done")



#=
@time p_DEFlux = DiffEqFlux.sciml_train(loss, p_initial, Fminbox(LBFGS()),
										lower_bounds=lower_bound,
										upper_bounds=upper_bound,
										show_trace=true,
										iterations=10,
										outer_iterations=10,
										f_tol=1e-4,
										outer_f_tol=5e-4,
										outer_x_tol=1e-6)
=#
println("True p: ",p_true)
println("Initial p: ", p_initial," Initial error: ",p_initial-p_true)
println("\nOptim Fminbox-LBFGS: ",sol_optim_minbox_lbfgs.u," \nLoss:",sol_optim_minbox_lbfgs.minimum," \nError: ",sol_optim_minbox_lbfgs.u-p_true)
println("\nNLopt Nelder-Mead: ",sol_nlopt_NedlerMead.u," \nLoss:",sol_nlopt_NedlerMead.minimum," \nError: ",sol_nlopt_NedlerMead.u-p_true)
println("\nNLopt LBFGS with Zygote gradient: ",sol_nlopt_LBFGS.u," \nLoss:",sol_nlopt_LBFGS.minimum," \nError: ",sol_nlopt_LBFGS.u-p_true)
println("\nNLopt LBFGS with FD gradient: ",sol_nlopt_LBFGS_fd.u," \nLoss:",sol_nlopt_LBFGS_fd.minimum," \nError: ",sol_nlopt_LBFGS_fd.u-p_true)
#println("\nNLopt LBFGS with ReverseDiff gradient: ",sol_nlopt_LBFGS_rd.u," \nLoss:",sol_nlopt_LBFGS_rd.minimum," \nError: ",sol_nlopt_LBFGS_rd.u-p_true)




#######################################################################

rho_Ramsey_u_optim,rho_Ramsey_v_optim = GLOQ.RamseyForwardSolve(rho_u0,rho_v0,
				 (2*pi).*sol_nlopt_LBFGS.u[1:3],omr,
				 [sol_nlopt_LBFGS.u[4:5];0.0],[sol_nlopt_LBFGS.u[6:7];0.0],
				 1,# initial state
				 TC,t_dark_times,N_states;initial_type = "states")
population_optim = GLOQ.get_population(rho_Ramsey_u_optim)


plot!(fig2,t_dark_times./1000.0,population_optim-population_Ramsey,
	  line=(:dash),lab=:false,legend=:outerright,)
display(fig2)

fig3=plot(t_dark_times./1000.0,population_optim-population_Ramsey,
	  line=(:dash),lab=:false,legend=:outerright,)
display(fig3)

fig4=plot(t_dark_times./1000.0,population_Ramsey)
plot!(fig4,t_dark_times./1000.0,population_optim,line=(:dash))
