"""
    RotationFrameDiagonal(omega::Array{Float64},omega_drive::Float64)

# Argument:
- omega: transition frequencies
- omega_drive: driving frequency
- N: number of states

# Output:
- Diagonal matrix diag{0,omega[1]-omr,omega[2]-omr,...,omega[N-1]-omr}
"""
function RotationFrameDiagonal(omega::Array{Float64},omega_drive::Float64)
	diag_vec = 0.0
	for i = 1:length(omega)
		# potential troublemaker of auto-diff
		diag_vec = [diag_vec;diag_vec[end]+omega[i]-omega_drive]
	end
	HK = Diagonal(diag_vec)
    return HK
end

function RotationFrameDiagonal(omega,omega_drive)
	diag_vec = 0.0
	for i = 1:length(omega)
		# potential troublemaker of auto-diff
		diag_vec = [diag_vec;diag_vec[end]+omega[i]-omega_drive]
	end
	HK = Diagonal(diag_vec)
    return HK
end

"""
    RotationFrameDiagonal(gamma1::Array{Float64},gamma2::Array{Float64})

# Argument:
- gamma1: determine the decay operator in the Lindblad system
- gamma2: determine the dephasing operator in the Lindblad system

# Output:
- L_decay: L_decay[i][i+1] = sqrt(gamma1[i]), and zeros everywhere else
- L_dephase: Diagonal{0,sqrt(gamma2[1]),...)
"""
function RotationFrameLindblad(gamma1::Array{Float64},gamma2::Array{Float64})
	#L_decay = Bidiagonal(zeros(length(gamma1)+1),sqrt.(gamma1),:U)
	L_dephase = Diagonal([0;sqrt.(gamma2)])
	L_decay = [zeros(Float64,length(gamma1)) Diagonal(sqrt.(gamma1))]
	L_decay = [L_decay;zeros(Float64,1,length(gamma1)+1)]
	#L_decay = [0.0 sqrt(gamma1[1]) 0.0 0.0;
	#		   0.0 0.0      sqrt(gamma1[2]) 0.0;
	#		   0.0 0.0      0.0      sqrt(gamma1[3]);
	#		   0.0 0.0      0.0      0.0]
    return L_decay, L_dephase
end

function RotationFrameLindblad(gamma1,gamma2)
	#L_decay = Bidiagonal(zeros(length(gamma1)+1),sqrt.(gamma1),:U)
	L_dephase = Diagonal([0;sqrt.(gamma2)])
	L_decay = [zeros(Float64,length(gamma1)) Diagonal(sqrt.(gamma1))]
	L_decay = [L_decay;zeros(Float64,1,length(gamma1)+1)]
	#L_decay = [0.0 sqrt(gamma1[1]) 0.0 0.0;
	#		   0.0 0.0      sqrt(gamma1[2]) 0.0;
	#		   0.0 0.0      0.0      sqrt(gamma1[3]);
	#		   0.0 0.0      0.0      0.0]
    return L_decay, L_dephase
end

"""
    RotationFrameRamseyControl(abs_control::Float64,phase_control::Float64,N::Int64)

# Argument:
- abs_control: strength of the control ``\\Omega``
- phase_control: phase of the control signal ``\\theta``
- N: number of states

# Output:
- HCK,HCS: the real and imaginary part of the control operator, determined by
  ``p(a+a^\\dagger)+iq(a-a^\\dagger),``
  where ``a`` is the lowering operator, ``p=|\\Omega|\\cos(\\theta)`` and ``q=|\\Omega|\\sin(\\theta)``
"""
function RotationFrameRamseyControl(abs_control::Float64,phase_control::Float64,N::Int64)
    p = abs_control*cos(phase_control);
    q = abs_control*sin(phase_control);
	amat = [Matrix(Diagonal(sqrt.(collect(1:N-1))));zeros(Float64,1,N-1)]
	#amat = [Diagonal(sqrt.(collect(1:N-1)));zeros(Float64,1,N-1)]
	amat = [zeros(Float64,N) amat]
    #amat = Bidiagonal(zeros(N),sqrt.(collect(1:N-1)),:U) # trouble maker, maybe we can make a Bidiagonal adjoint
    HCK = p*(amat+transpose(amat))
    HCS = q*(amat-transpose(amat))
    return HCK,HCS
end
"""
    RotationFrameRamseyControl(N::Int64)

# Argument:
- abs_control: strength of the control ``\\Omega``
- phase_control: phase of the control signal ``\\theta``
- N: number of states

# Output:
- HCK,HCS: the real and imaginary part of the control operator, determined by
  ``(a+a^\\dagger)``
"""
function RotationFrameRamseyControl(N::Int64)
	amat = [Matrix(Diagonal(sqrt.(collect(1:N-1))));zeros(Float64,1,N-1)]
	amat = [zeros(Float64,N) amat]
    HCK = amat+transpose(amat)
	HCS = zeros(Float64,N,N)
    return HCK,HCS
end

"""
function RamseyForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
						  omega::Array{Float64},omega_drive::Float64,
						  gamma1::Array{Float64},gamma2::Array{Float64},
						  InitialState::Int64,
						  TC::Float64,t_dark_times::Array{Float64},
						  N_states::Int64=0;
						  initial_type="states",
						  method="exponential",
						  DiffEqKwargs...)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states
- InitialState: initial state of the density matrix
- method: method to solve the Lindblad system
 1. "exponential": exponential time integrator
 2. "DiffEqDefault": a default choice made by DifferentialEquations.jl
 3. Other solvers availabe in DifferentialEquations.jl, for example, method = Trapezoid()
- DiffEqKwargs: keyword arguments feed to the ``solve function" of DifferentialEquations.jl

# Output:
- rho_ramsey_u,rho_ramsey_v: density matrix at dark times, with ``\\rho=\\rho_u-i\\rho_v``
"""
function RamseyForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   		  omega::Array{Float64},omega_drive::Float64,
				   		  gamma1::Array{Float64},gamma2::Array{Float64},
						  InitialState::Int64,
				   		  TC::Float64,t_dark_times::Array{Float64},
				   		  N_states::Int64=0;
						  initial_type="states",
						  method="exponential",
						  DiffEqKwargs...)
    if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	if(method=="exponential")
		Ω = 0.5*pi/(TC*sqrt(InitialState+1.0))
		θ = 0.0
		HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
	else
		HK_control,HS_control = RotationFrameRamseyControl(N_states)
	end
	# Lindblad operator
	L1,L2 = RotationFrameLindblad(gamma1,gamma2)
	# Assemble the operators for the vectorized system
	# Half pi pulse
	LK_half_pi,LS_half_pi,LD = make_lindblad_operator(Matrix(HK_control+HK_free),Matrix(HS_control),(L1,L2))
	# Free evolution
	LK_free = make_hamiltonian_operator(Matrix(HK_free))
	#################################
	# Step 2: Forward solves
	#################################
	# initial conditions
	if(initial_type == "states")
		rho_u_initial = (rho_u0*transpose(rho_u0)+rho_v0*transpose(rho_v0))[:]
		rho_v_initial = (rho_v0*transpose(rho_u0)-rho_u0*transpose(rho_v0))[:]
		rho_vec0 = [rho_u_initial;rho_v_initial]
	elseif(initial_type == "density")
		rho_vec0 = [rho_u0;rho_v0];
	else
		println("Error! initial_type must be \"density\" or \"states\"")
		return
	end
	N_dark_times = length(t_dark_times)
	if(method=="exponential")
	#if(false)
		# Ramsey experiment
		Half_pi_operator = exp(0.5*TC*[LD+LS_half_pi -LK_half_pi; LK_half_pi LD+LS_half_pi])
		FreeEvolution = [ LD -LK_free;
					  	  LK_free  LD]
		# half-pi pusle
		rho_vec1 = Half_pi_operator*rho_vec0
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# half-pi pulse
		rho_ramsey = Half_pi_operator*rho_vec2
		for i = 2:N_dark_times
			# free propagation for a dark time
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_vec2 = Free_operator*rho_vec1
			# half-pi pulse
			rho_ramsey = [rho_ramsey Half_pi_operator*rho_vec2]
		end
	else
		function GaussianPulse(_time,_amplitude,_width)
			σ = _width/sqrt(2.0*pi)
			return _amplitude*exp(-(_time-0.5*TC)^2/(2.0*σ*σ))
		end
		width = TC/2.5
		GaussianAmplitude = 0.5*pi/(width*sqrt(InitialState+1.0))
		function ControlledODE!(_drho,_rho,_p,_t)
			# Half pi pulse
			_LK_half_pi,_LS_half_pi,_LD = make_lindblad_operator(HK_control*GaussianPulse(_t,GaussianAmplitude,width)+HK_free,
															     HS_control*GaussianPulse(_t,GaussianAmplitude,width),
															     (L1,L2))
			SystemOperator  = [_LD+_LS_half_pi -_LK_half_pi;
					  		   _LK_half_pi      _LD+_LS_half_pi]
			_drho .= SystemOperator*_rho
		end
		# DiffEq interface for the controlled problem
		controlled_problem = ODEProblem(ControlledODE!,rho_vec0,(0.0,0.5*TC),[1.0])
		# Free propagation part
		FreeEvolution = [ LD -LK_free;
						  LK_free  LD]
		# half-pi pusle
		if(method=="DiffEqDefault")
			sol = solve(controlled_problem;DiffEqKwargs...)
		else
			sol = solve(controlled_problem,method;DiffEqKwargs...)
		end
		rho_vec1 = sol[:,end]
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# second half-pi pulse
		controlled_problem = remake(controlled_problem,u0=rho_vec2)
		if(method=="DiffEqDefault")
			sol = solve(controlled_problem;DiffEqKwargs...)
		else
			sol = solve(controlled_problem,method;DiffEqKwargs...)
		end
		rho_ramsey = sol[:,end]
		# free propagation for a dark time + second half-pi pulse
		n_dark_times = length(t_dark_times)
		for i = 2:n_dark_times
			# free propagation
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_vec2 = Free_operator*rho_vec1
			# second half pi pulse
			controlled_problem = remake(controlled_problem,u0=rho_vec2)
			if(method=="DiffEqDefault")
				sol = solve(controlled_problem;DiffEqKwargs...)
			else
				sol = solve(controlled_problem,method;DiffEqKwargs...)
			end
			rho_ramsey = [ rho_ramsey sol[:,end] ]
		end
	end
	return rho_ramsey[1:N_states^2,:],rho_ramsey[N_states^2+1:end,:]
end

# Defined for the compatability of adjoint solves in Turing.jl to hold their data structure
function RamseyForwardSolve(rho_u0,rho_v0,
				   		  omega,omega_drive,
				   		  gamma1,gamma2,
						  InitialState,
				   		  TC,t_dark_times,
				   		  N_states::Int64=0;
						  initial_type="states",
						  method="exponential",
						  DiffEqKwargs...)
    if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	if(method=="exponential")
		Ω = 0.5*pi/(TC*sqrt(InitialState+1.0))
		θ = 0.0
		HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
	else
		HK_control,HS_control = RotationFrameRamseyControl(N_states)
	end
	# Lindblad operator
	L1,L2 = RotationFrameLindblad(gamma1,gamma2)
	# Assemble the operators for the vectorized system
	# Half pi pulse
	LK_half_pi,LS_half_pi,LD = make_lindblad_operator(Matrix(HK_control+HK_free),Matrix(HS_control),(L1,L2))
	# Free evolution
	LK_free = make_hamiltonian_operator(Matrix(HK_free))
	#################################
	# Step 2: Forward solves
	#################################
	# initial conditions
	if(initial_type == "states")
		rho_u_initial = (rho_u0*transpose(rho_u0)+rho_v0*transpose(rho_v0))[:]
		rho_v_initial = (rho_v0*transpose(rho_u0)-rho_u0*transpose(rho_v0))[:]
		rho_vec0 = [rho_u_initial;rho_v_initial]
	elseif(initial_type == "density")
		rho_vec0 = [rho_u0;rho_v0];
	else
		println("Error! initial_type must be \"density\" or \"states\"")
		return
	end
	N_dark_times = length(t_dark_times)
	if(method=="exponential")
	#if(false)
		# Ramsey experiment
		Half_pi_operator = exp(0.5*TC*[LD+LS_half_pi -LK_half_pi; LK_half_pi LD+LS_half_pi])
		FreeEvolution = [ LD -LK_free;
					  	  LK_free  LD]
		# half-pi pusle
		rho_vec1 = Half_pi_operator*rho_vec0
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# half-pi pulse
		rho_ramsey = Half_pi_operator*rho_vec2
		for i = 2:N_dark_times
			# free propagation for a dark time
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_vec2 = Free_operator*rho_vec1
			# half-pi pulse
			rho_ramsey = [rho_ramsey Half_pi_operator*rho_vec2]
		end
	else
		function GaussianPulse(_time,_amplitude,_width)
			σ = _width/sqrt(2.0*pi)
			return _amplitude*exp(-(_time-0.5*TC)^2/(2.0*σ*σ))
		end
		width = TC/2.5
		GaussianAmplitude = 0.5*pi/(width*sqrt(InitialState+1.0))
		function ControlledODE!(_drho,_rho,_p,_t)
			# Half pi pulse
			_LK_half_pi,_LS_half_pi,_LD = make_lindblad_operator(HK_control*GaussianPulse(_t,GaussianAmplitude,width)+HK_free,
															     HS_control*GaussianPulse(_t,GaussianAmplitude,width),
															     (L1,L2))
			SystemOperator  = [_LD+_LS_half_pi -_LK_half_pi;
					  		   _LK_half_pi      _LD+_LS_half_pi]
			_drho .= SystemOperator*_rho
		end
		# DiffEq interface for the controlled problem
		controlled_problem = ODEProblem(ControlledODE!,rho_vec0,(0.0,0.5*TC),[omega;gamma1;gamma2])
		# Free propagation part
		FreeEvolution = [ LD -LK_free;
						  LK_free  LD]
		# half-pi pusle
		if(method=="DiffEqDefault")
			sol = solve(controlled_problem;DiffEqKwargs...)
		else
			sol = solve(controlled_problem,method;DiffEqKwargs...)
		end
		rho_vec1 = sol[:,end]
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# second half-pi pulse
		controlled_problem = remake(controlled_problem,u0=rho_vec2)
		if(method=="DiffEqDefault")
			sol = solve(controlled_problem;DiffEqKwargs...)
		else
			sol = solve(controlled_problem,method;DiffEqKwargs...)
		end
		rho_ramsey = sol[:,end]
		# free propagation for a dark time + second half-pi pulse
		n_dark_times = length(t_dark_times)
		for i = 2:n_dark_times
			# free propagation
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_vec2 = Free_operator*rho_vec1
			# second half pi pulse
			controlled_problem = remake(controlled_problem,u0=rho_vec2)
			if(method=="DiffEqDefault")
				sol = solve(controlled_problem;DiffEqKwargs...)
			else
				sol = solve(controlled_problem,method;DiffEqKwargs...)
			end
			rho_ramsey = [ rho_ramsey sol[:,end] ]
		end
	end
	return rho_ramsey[1:N_states^2,:],rho_ramsey[N_states^2+1:end,:]
end
"""
    RamseyParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive:Float64,
					   charge_noise::Array{Float64},
					   gamma1::Array{Float64},gamma2::Array{Float64},
					   InitalState::Int64,
					   TC::Float64,t_dark_times::Array{Float64},
					   N_states::Int64=0;initial_type="states")

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- charge_noise: charge noise
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- InitialState: initial state of the density matrix
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states
- InitialState: initial state of the density matrix
- method: method to solve the Lindblad system
 1. "exponential": exponential time integrator
 2. "DiffEqDefault": a default choice made by DifferentialEquations.jl
 3. Other solvers availabe in DifferentialEquations.jl, for example, method = Trapezoid()
- DiffEqKwargs: keyword arguments feed to the ``solve function" of DifferentialEquations.jl

# Parity event
- Transition frequency will take equal probability to be omega``\\pm``0.5``\\times``charge_noise

# Output:
- rho_ramsey_u,rho_ramsey_v: density matrix at dark times, with ``\\rho=\\rho_u-i\\rho_v``
"""
function RamseyParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   		  omega::Array{Float64},omega_drive::Float64,
						  charge_noise::Array{Float64},
				   		  gamma1::Array{Float64},gamma2::Array{Float64},
						  InitialState::Int64,
				   		  TC::Float64,t_dark_times::Array{Float64},
				   		  N_states::Int64=0;
						  initial_type="states",
						  method="exponential",
						  DiffEqKwargs...)

	rho_u1, rho_v1 = RamseyForwardSolve(rho_u0,rho_v0,
					   		  omega-0.5.*charge_noise,omega_drive,
					   		  gamma1,gamma2,
							  InitialState,
					   		  TC,t_dark_times,
					   		  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	rho_u2, rho_v2 = RamseyForwardSolve(rho_u0,rho_v0,
							  omega+0.5.*charge_noise,omega_drive,
							  gamma1,gamma2,
							  InitialState,
							  TC,t_dark_times,
							  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	return 0.5.*(rho_u1+rho_u2),0.5.*(rho_v1+rho_v2)
end

"""
    EchoForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive::Float64,
					   gamma1::Array{Float64},gamma2::Array{Float64},
					   InitialState::Int64,
					   TC::Float64,t_dark_times::Array{Float64},
					   N_states::Int64=0;
					   initial_type="states",
					   method="exponential",
					   DiffEqKwargs...)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states
- InitialState: initial state of the density matrix
- method: method to solve the Lindblad system
 1. "exponential": exponential time integrator
 2. "DiffEqDefault": a default choice made by DifferentialEquations.jl
 3. Other solvers availabe in DifferentialEquations.jl, for example, method = Trapezoid()
- DiffEqKwargs: keyword arguments feed to the ``solve function" of DifferentialEquations.jl
"""
function EchoForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   omega::Array{Float64},omega_drive::Float64,
				   gamma1::Array{Float64},gamma2::Array{Float64},
				   InitialState::Int64,
				   TC::Float64,t_dark_times::Array{Float64},
				   N_states::Int64=0;
				   initial_type="states",
				   method="exponential",
				   DiffEqKwargs...)

	if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	if(method=="exponential")
		Ω = 0.5*pi/(TC*sqrt(InitialState+1.0))
		θ = 0.0
		HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
	else
		HK_control,HS_control = RotationFrameRamseyControl(N_states)
	end
	# Lindblad operator
	L1,L2 = RotationFrameLindblad(gamma1,gamma2)
	# Assemble the operators for the vectorized system
	# Half pi pulse
	LK_half_pi,LS_half_pi,LD = make_lindblad_operator(Matrix(HK_control+HK_free),Matrix(HS_control),(L1,L2))
	# Free evolution
	LK_free = make_hamiltonian_operator(Matrix(HK_free))
	#################################
	# Step 2: Forward solves
	#################################
	# initial conditions
	if(initial_type == "states")
		rho_u_initial = (rho_u0*transpose(rho_u0)+rho_v0*transpose(rho_v0))[:]
		rho_v_initial = (rho_v0*transpose(rho_u0)-rho_u0*transpose(rho_v0))[:]
		rho_vec0 = [rho_u_initial;rho_v_initial]
	elseif(initial_type == "density")
		rho_vec0 = [rho_u0;rho_v0];
	else
		println("Error! initial_type must be \"density\" or \"states\"")
		return
	end
	N_dark_times = length(t_dark_times)
	if(method=="exponential")
		# Assemble operators
		Half_pi_operator = exp(0.5*TC*[LD+LS_half_pi -LK_half_pi; LK_half_pi LD+LS_half_pi])
		Pi_operator = Half_pi_operator*Half_pi_operator
		FreeEvolution = [ LD -LK_free;
					      LK_free  LD]
					  # half-pi pusle
		rho_vec1 = Half_pi_operator*rho_vec0
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# pi pulse
		rho_vec2 = Pi_operator*rho_vec2
		# free propagation for the same dark time
		rho_vec2 = Free_operator*rho_vec2
		# half-pi pulse
		rho_echo = Half_pi_operator*rho_vec2
		for i = 2:N_dark_times
			# free propagation for a dark time
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_vec2 = Free_operator*rho_vec1
			# pi pulse
			rho_vec2 = Pi_operator*rho_vec2
			# free propagation for the same dark time
			rho_vec2 = Free_operator*rho_vec2
			# half-pi pulse
			rho_echo = [rho_echo Half_pi_operator*rho_vec2]
		end
	else
		function GaussianPulse(_time,_amplitude,_width)
			_tot_time = _width*2.5
			σ = _width/sqrt(2.0*pi)
			return _amplitude*exp(-(_time-0.5*_tot_time)^2/(2.0*σ*σ))
		end
		width = TC/2.5
		GaussianAmplitude = 0.5*pi/(width*sqrt(InitialState+1.0))
		function ControlledODE!(_drho,_rho,_p,_t)
			# Half pi pulse
			_LK_half_pi,_LS_half_pi,_LD = make_lindblad_operator(HK_control*GaussianPulse(_t,GaussianAmplitude,width)+HK_free,
															     HS_control*GaussianPulse(_t,GaussianAmplitude,width),
															     (L1,L2))
			SystemOperator  = [_LD+_LS_half_pi -_LK_half_pi;
					  		   _LK_half_pi      _LD+_LS_half_pi]
			_drho .= SystemOperator*_rho
		end
		# DiffEq interface for the controlled problem
		half_pi_pulse = ODEProblem(ControlledODE!,rho_vec0,(0.0,0.5*TC),[1.0])
		# Free propagation part
		FreeEvolution = [ LD -LK_free;
						  LK_free  LD]
		# half-pi pusle
		if(method=="DiffEqDefault")
			sol = solve(half_pi_pulse;DiffEqKwargs...)
		else
			sol = solve(half_pi_pulse,method;DiffEqKwargs...)
		end
		rho_vec1 = sol[:,end]
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_vec2 = Free_operator*rho_vec1
		# pi pulse
		pi_pulse = ODEProblem(ControlledODE!,rho_vec2,(0.0,TC),[1.0])
		if(method=="DiffEqDefault")
			sol = solve(pi_pulse;DiffEqKwargs...)
		else
			sol = solve(pi_pulse,method;DiffEqKwargs...)
		end
		# free propagation for the same dark time
		rho_vec2 = Free_operator*sol[:,end]
		# half pi pulse
		half_pi_pulse = remake(half_pi_pulse,u0=rho_vec2)
		if(method=="DiffEqDefault")
			sol = solve(half_pi_pulse;DiffEqKwargs...)
		else
			sol = solve(half_pi_pulse,method;DiffEqKwargs...)
		end
		rho_echo = sol[:,end]
		# free propagation for a dark time + second half-pi pulse
		n_dark_times = length(t_dark_times)
		for i = 2:n_dark_times
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			# second half pi
			rho_vec2 = Free_operator*rho_vec1
			# pi pulse
			pi_pulse = remake(pi_pulse,u0=rho_vec2)
			if(method=="DiffEqDefault")
				sol = solve(pi_pulse;DiffEqKwargs...)
			else
				sol = solve(pi_pulse,method;DiffEqKwargs...)
			end
			# free propagation for the same dark time
			rho_vec2 = Free_operator*sol[:,end]
			# half pi pulse
			half_pi_pulse = remake(half_pi_pulse,u0=rho_vec2)
			if(method=="DiffEqDefault")
				sol = solve(half_pi_pulse;DiffEqKwargs...)
			else
				sol = solve(half_pi_pulse,method;DiffEqKwargs...)
			end
			rho_echo = [rho_echo sol[:,end]]
		end

	end

	return rho_echo[1:N_states^2,:],rho_echo[N_states^2+1:end,:]
end

"""
    EchoParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive::Float64,
					   charge_noise::Array{Float64},
					   gamma1::Array{Float64},gamma2::Array{Float64},
					   InitialState::Int64,
					   TC::Float64,t_dark_times::Array{Float64},
					   N_states::Int64=0;
					   initial_type="states",
					   method="exponential",
					   DiffEqKwargs...)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- charge_noise: charge noise
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- InitialState: initial state of the density matrix
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states
- InitialState: initial state of the density matrix
- method: method to solve the Lindblad system
	- "exponential": exponential time integrator
	- "DiffEqDefault": a default choice made by DifferentialEquations.jl
	- Other solvers availabe in DifferentialEquations.jl, for example, method = Trapezoid()
- DiffEqKwargs: keyword arguments feed to the ``solve function" of DifferentialEquations.jl

# Parity event
- Transition frequency will take equal probability to be omega``\\pm``0.5``\\times``charge_noise

# Output:
- rho_echo_u,rho_echo_v: density matrix at dark times, with ``\\rho=\\rho_u-i\\rho_v``
"""
function EchoParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   omega::Array{Float64},omega_drive::Float64,
				   charge_noise::Array{Float64},
				   gamma1::Array{Float64},gamma2::Array{Float64},
				   InitialState::Int64,
				   TC::Float64,t_dark_times::Array{Float64},
				   N_states::Int64=0;
				   initial_type="states",
				   method="exponential",
				   DiffEqKwargs...)

	rho_u1, rho_v1 = EchoForwardSolve(rho_u0,rho_v0,
							  omega-0.5.*charge_noise,omega_drive,
							  gamma1,gamma2,
							  InitialState,
							  TC,t_dark_times,
							  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	rho_u2, rho_v2 = EchoForwardSolve(rho_u0,rho_v0,
							  omega+0.5.*charge_noise,omega_drive,
							  gamma1,gamma2,
							  InitialState,
							  TC,t_dark_times,
							  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	return 0.5.*(rho_u1+rho_u2),0.5.*(rho_v1+rho_v2)
end

"""
    T1ForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive::Float64,
					   gamma1::Array{Float64},gamma2::Array{Float64},
					   InitialState::Int64,
					   TC::Float64,t_dark_times::Array{Float64},
					   N_states::Int64=0;
					   initial_type="states",
					   method="exponential",
					   DiffEqKwargs...)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- charge_noise: charge noise
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- InitialState: initial state of the density matrix
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states

# Output:
- rho_echo_u,rho_echo_v: density matrix at dark times, with ``\\rho=\\rho_u-i\\rho_v``
"""
function T1ForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   omega::Array{Float64},omega_drive::Float64,
				   gamma1::Array{Float64},gamma2::Array{Float64},
				   InitialState::Int64,
				   TC::Float64,t_dark_times::Array{Float64},
				   N_states::Int64=0;
				   initial_type="states",
				   method="exponential",
				   DiffEqKwargs...)

	if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	# Hamiltonain with control
	if(method=="exponential")
		Ω = 0.5*pi/(TC*sqrt(InitialState+1.0))
		θ = 0.0
		HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
	else
		HK_control,HS_control = RotationFrameRamseyControl(N_states)
	end
	# Lindblad operator
	L1,L2 = RotationFrameLindblad(gamma1,gamma2)
	# Assemble the operators for the vectorized system
	# Half pi pulse
	LK_half_pi,LS_half_pi,LD = make_lindblad_operator(Matrix(HK_control+HK_free),Matrix(HS_control),(L1,L2))
	# Free evolution
	LK_free = make_hamiltonian_operator(Matrix(HK_free))
	#################################
	# Step 2: Forward solves
	#################################
	# initial conditions
	if(initial_type == "states")
		rho_u_initial = (rho_u0*transpose(rho_u0)+rho_v0*transpose(rho_v0))[:]
		rho_v_initial = (rho_v0*transpose(rho_u0)-rho_u0*transpose(rho_v0))[:]
		rho_vec0 = [rho_u_initial;rho_v_initial]
	elseif(initial_type == "density")
		rho_vec0 = [rho_u0;rho_v0];
	else
		println("Error! initial_type must be \"density\" or \"states\"")
		return
	end
	N_dark_times = length(t_dark_times)
	if(method=="exponential")
		# Assemble operators
		Pi_operator = exp(TC*[LD+LS_half_pi -LK_half_pi; LK_half_pi LD+LS_half_pi])
		FreeEvolution = [ LD -LK_free;
					  	  LK_free  LD]
		# pi pusle
		rho_vec1 = Pi_operator*rho_vec0
		# free propagation for a dark time
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_T1 = Free_operator*rho_vec1
		for i = 2:N_dark_times
			# free propagation for a dark time
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_T1 = [rho_T1 Free_operator*rho_vec1]
		end
	else
		function GaussianPulse(_time,_amplitude,_width)
			σ = _width/sqrt(2.0*pi)
			return _amplitude*exp(-(_time-0.5*TC)^2/(2.0*σ*σ))
		end
		width = TC/2.5
		GaussianAmplitude = 0.5*pi/(width*sqrt(InitialState+1.0))
		function ControlledODE!(_drho,_rho,_p,_t)
			# Half pi pulse
			_LK_half_pi,_LS_half_pi,_LD = make_lindblad_operator(HK_control*GaussianPulse(_t,GaussianAmplitude,width)+HK_free,
																 HS_control*GaussianPulse(_t,GaussianAmplitude,width),
																 (L1,L2))
			SystemOperator  = [_LD+_LS_half_pi -_LK_half_pi;
							   _LK_half_pi      _LD+_LS_half_pi]
			_drho .= SystemOperator*_rho
		end
		# DiffEq interface for the controlled problem
		pi_pulse = ODEProblem(ControlledODE!,rho_vec0,(0.0,TC),[1.0])
		# Free propagation part
		FreeEvolution = [ LD -LK_free;
						  LK_free  LD]
		# pi pusle
		if(method=="DiffEqDefault")
			sol = solve(pi_pulse;DiffEqKwargs...)
		else
			sol = solve(pi_pulse,method;DiffEqKwargs...)
		end
		rho_vec1 = sol[:,end]
		Free_operator = exp(t_dark_times[1]*FreeEvolution)
		rho_T1 = Free_operator*rho_vec1
		for i = 2:N_dark_times
			# free propagation for a dark time
			Free_operator = exp(t_dark_times[i]*FreeEvolution)
			rho_T1 = [rho_T1 Free_operator*rho_vec1]
		end
	end
	return rho_T1[1:N_states^2,:],rho_T1[N_states^2+1:end,:]

end

"""
    T1ParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive::Float64,
					   charge_noise::Array{Float64},
					   gamma1::Array{Float64},gamma2::Array{Float64},
					   InitialState::Int64,
					   TC::Float64,t_dark_times::Array{Float64},
					   N_states::Int64=0;
					   initial_type="states",
					   method="exponential",
					   DiffEqKwargs...)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- charge_noise: charge noise
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- InitialState: initial state of the density matrix
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states
- InitialState: initial state of the density matrix
- method: method to solve the Lindblad system
	- "exponential": exponential time integrator
	- "DiffEqDefault": a default choice made by DifferentialEquations.jl
	- Other solvers availabe in DifferentialEquations.jl, for example, method = Trapezoid()
- DiffEqKwargs: keyword arguments feed to the ``solve function" of DifferentialEquations.jl

# Parity event
- Transition frequency will take equal probability to be omega``\\pm``0.5``\\times``charge_noise

# Output:
- rho_echo_u,rho_echo_v: density matrix at dark times, with ``\\rho=\\rho_u-i\\rho_v``
"""
function T1ParityForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   omega::Array{Float64},omega_drive::Float64,
				   charge_noise::Array{Float64},
				   gamma1::Array{Float64},gamma2::Array{Float64},
				   InitialState::Int64,
				   TC::Float64,t_dark_times::Array{Float64},
				   N_states::Int64=0;
				   initial_type="states",
				   method="exponential",
				   DiffEqKwargs...)

	rho_u1, rho_v1 = T1ForwardSolve(rho_u0,rho_v0,
							  omega-0.5.*charge_noise,omega_drive,
							  gamma1,gamma2,
							  InitialState,
							  TC,t_dark_times,
							  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	rho_u2, rho_v2 = T1ForwardSolve(rho_u0,rho_v0,
							  omega+0.5.*charge_noise,omega_drive,
							  gamma1,gamma2,
							  InitialState,
							  TC,t_dark_times,
							  N_states;
							  initial_type,
							  method,
							  DiffEqKwargs...)
	return 0.5.*(rho_u1+rho_u2),0.5.*(rho_v1+rho_v2)
end
