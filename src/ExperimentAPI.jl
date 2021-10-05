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
	omega_sum = 0.0
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
	amat = [Diagonal(sqrt.(collect(1:N-1)));zeros(Float64,1,N-1)]
	amat = [zeros(Float64,N) amat]
    #amat = Bidiagonal(zeros(N),sqrt.(collect(1:N-1)),:U) # trouble maker, maybe we can make a Bidiagonal adjoint
    HCK = p*(amat+transpose(amat))
    HCS = q*(amat-transpose(amat))
    return HCK,HCS
end

"""
    RamseyForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
					   omega::Array{Float64},omega_drive:Float64,
					   gamma1::Array{Float64},gamma2::Array{Float64}
					   TC::Float64,t_dark_times::Float64,
					   N_states::Int64=0)

# Argument:
- rho_u0,rho_v0: initial states ``\\rho_{u_0}-i\\rho_{v_0}``
- omega: transition frequencies
- omega_drive: driving frequency
- gamma1,gamma2: determine the decay part and the dephasing part of the Lindblad operators
- TC: control time of the control signal
- t_dark_times: the dark time in the Ramsey experiment
- N_states: number of states

# Output:
-
"""
function RamseyForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   		  omega::Array{Float64},omega_drive::Float64,
				   		  gamma1::Array{Float64},gamma2::Array{Float64},
				   		  TC::Float64,t_dark_times::Array{Float64},
				   		  N_states::Int64=0;initial_type="states")
    if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	Ω = 0.5*pi/(TC*sqrt(2.0))
	θ = 0.0
	HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
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

	return rho_ramsey[1:N_states^2,:],rho_ramsey[N_states^2+1:end,:]
end



function RamseyForwardSolve(rho_u0::Array{Float64},rho_v0::Array{Float64},
				   		  omega::Array{Float64},omega_drive::Float64,
				   		  gamma1::Array{Float64},gamma2::Array{Float64},
				   		  TC::Float64,t_dark_times::Array{Float64},
				   		  N_states::Int64=0;initial_type="states")
    if(N_states==0)
		N_states = length(omega)+1
	end
	#################################
	# Step 1: assemble operators
	#################################
	# Hamiltonian without control
	HK_free = RotationFrameDiagonal(omega,omega_drive)
	# Hamiltonain with control
	Ω = 0.5*pi/(TC*sqrt(2.0))
	θ = 0.0
	HK_control,HS_control = RotationFrameRamseyControl(Ω,θ,N_states)
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

	return rho_ramsey[1:N_states^2,:],rho_ramsey[N_states^2+1:end,:]
end
