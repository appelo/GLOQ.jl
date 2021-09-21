using LinearAlgebra
using DifferentialEquations
##################################################
# Exponential solver
# input: L, the whole propagate operator, t_span
# where the function will be evaluated, u0
# initial states
# output: corresponding solutions
##################################################
function exponential_solver(rho_vec0,L,t_span::Array{Float64};initial_type = "density")
    if(initial_type == "states")
        rho_vec0 = (rho_vec0*rho_vec0')[:]
    elseif(initial_type != "density")
        println("Error! initial_type must be \"density\" or \"states\"")
        return
    end
    rho_vec = zeros(ComplexF64,length(rho_vec0),length(t_span))
    for time_step = 1:length(t_span)
        rho_vec[:,time_step] = exp(L*t_span[time_step])*rho_vec0
    end
    return rho_vec
end

function exponential_solver(rho_vec0,L,
                            t_span::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}};
                            initial_type = "density")
    if(initial_type == "states")
        rho_vec0 = (rho_vec0*rho_vec0')[:]
    elseif(initial_type != "density")
        println("Error! initial_type must be \"density\" or \"states\"")
        return
    end
    rho_vec = zeros(ComplexF64,length(rho_vec0),length(t_span))
    dt = t_span[2]-t_span[1]
    Propagator = exp(L*dt)
    rho_vec[:,1] = exp(L*t_span[1])*rho_vec0
    for time_step = 2:length(t_span)
        rho_vec[:,time_step] = Propagator*rho_vec[:,time_step-1]
    end
    return rho_vec
end

##################################################
# Differential equations solver:
# provide interfaces to DifferentialEquations
# package.
# L: Lindblad operator
# time/tspan: final time or time beining evaluated
# rho0: initial condition
##################################################
# constant Lindblad problem
function LindbladODEProblem(rho0,L::Array{ComplexF64,2},time_final::Float64;initial_type = "density")
    if(initial_type=="states")
        rho0 = (rho0*rho0')[:];
    end
    rho0 = convert(Vector{ComplexF64},rho0)
    println(typeof(rho0))
    function ode_problem!(_drho,_rho,_p,_t)
    	# .= is essential, = will not work
    	_drho .= _p*_rho
    end
    return ODEProblem(ode_problem!,rho0,(0.0,time_final),L)
end
