"""
    exponential_solver_complex(rho_vec0,L,t_span::Array{Float64};initial_type = "density"):
# Purpose:
- use an exponential integrator to integrate

    ``\\rho_t=-i(L_K+iL_S) \\rho+L_D \rho``

# Argument:
- ``\\rho_{u0},\\rho_{v0}``: ``\\rho0 = ``\\rho_{u0} - i \\rho_{v0}``
- LK: real part of the Hamiltonian operator
- LS: imaginary part of the Hamiltonain operator
- LD: Lindblad operator
- t_span: where the solution will be stored
- initial_type: "density" vectorized density matrix, "states" states vector
# Output:
- rho_u,rho_v: solutions at t_span given by rho(t_i) = rho_u(:,i)-i rho_v(:,i)
"""
function exponential_solver(rho_u0,rho_v0,
                            LK,LS,LD,
                            t_span,N=0;
                            initial_type = "density")
    if(initial_type == "states")
        rho_u_initial = (rho_u0*transpose(rho_u0)+rho_v0*transpose(rho_v0))[:]
        rho_v_initial = (rho_v0*transpose(rho_u0)-rho_u0*transpose(rho_v0))[:]
        rho_vec0 = [rho_u_initial;rho_v_initial]
        if(N==0)
            N = length(rho_u_initial)
        end
    elseif(initial_type == "density")
        rho_vec0 = [rho_u0;rho_v0];
        if(N==0)
            N = length(rho_u0)
        end
    else
        println("Error! initial_type must be \"density\" or \"states\"")
        return
    end
    rho_vec = zeros(Float64,length(rho_vec0),length(t_span))
    LL = [LS+LD -LK;
          LK     LS+LD]
    for time_step = 1:length(t_span)
        rho_vec[:,time_step] = exp(LL*t_span[time_step])*rho_vec0
    end
    return rho_vec[1:N,:],rho_vec[N+1:2*N,:]
end
"""
    exponential_solver_complex(rho_vec0,L,t_span::Array{Float64};initial_type = "density"):

# Argument:
- L: the whole propagation operator
- t_span: where the function will be evaluated stored in
- rho0_vec: initial density/state
- initial_type: decide we aer given an initial density matrix or a state vector

# Output:
- solutions at t_span
"""
function exponential_solver_complex(rho_vec0,L,t_span::Array{Float64};initial_type = "density")
    if(initial_type == "states")
        rho_vec0 = (rho_vec0*adjoint(rho_vec0))[:]
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
"""
    exponential_solver(rho_vec0,L,t_span::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}};
                       initial_type = "density")

# Argument:
- L: the whole propagation operator
- t_span: where the function will be evaluated stored in
- rho0_vec: initial density/state
- initial_type: decide we aer given an initial density matrix or a state vector

# Output:
- solutions at t_span
"""
function exponential_solver_complex(rho_vec0,L,
                            t_span::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}};
                            initial_type = "density")
    if(initial_type == "states")
        rho_vec0 = (rho_vec0*adjoint(rho_vec0))[:]
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

"""
    LindbladODEProblemComplex(rho0,L::Array{ComplexF64,2},time_final::Float64;initial_type = "density")

Function provide interfaces to DifferentialEquations
package to solve: ``\\rho_t = L \\rho``

# Argument:
- L: Lindblad operator
- time: final time or time beining evaluated
- rho0: initial condition
- initial_type: specify initial value is a density matrix/a state vector

# Output:
A problem object which we will feed to DifferentialEquations.jl
"""
function LindbladODEProblemComplex(rho0,L::Array{ComplexF64,2},time_final::Float64;initial_type = "density")
    if(initial_type=="states")
        rho0 = (rho0*adjoint(rho0))[:];
    end
    rho0 = convert(Vector{ComplexF64},rho0)
    function ode_problem!(_drho,_rho,_p,_t)
    	_drho .= _p*_rho
    end
    return ODEProblem(ode_problem!,rho0,(0.0,time_final),L)
end

"""
    LindbladODEProblem(rho_u0,L::Array{ComplexF64,2},time_final::Float64;initial_type = "density")

Function provide interfaces to DifferentialEquations
package to solve the real-valued Lindblad system:
``(\\rho_u - i \\rho_v)_t = -i (L_K+iL_S)(\\rho_u-i \\rho_v) + L_D(\\rho_u-i \\rho_v)``

# Argument:
- ``\\rho_{u0}, \\rho_{v0}``: initial condition ``\\rho_0 = \\rho_{u0} - i \\rho_{v0}``
- LK: real part of the Hamiltonian operator
- LS: imaginary part of the Hamiltonain operator
- LD: Lindblad operator
- time: final time or time beining evaluated
- initial_type: specify initial value is a density matrix/a state vector

# Output:
A problem object which we will feed to DifferentialEquations.jl
"""
function LindbladODEProblem(rho_u0,rho_v0,LK::Array{Float64,2},LS::Array{Float64},LD::Array{Float64,2},time_final::Float64;initial_type = "density")
    if(initial_type=="states")
        rho_u0,rho_v0 = convert_state_to_density(rho_u0,rho_v0)
    end
    L_real = [LS+LD -LK;
              LK     LS+LD]
    function ode_problem!(_drho,_rho,_p,_t)
    	_drho .= _p*_rho
    end
    return ODEProblem(ode_problem!,[rho_u0;rho_v0],(0.0,time_final),L_real)
end
