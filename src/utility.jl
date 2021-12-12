"""
	get_population(rho_vec::Array{ComplexF64,1})

# Argument:
- rho_vec, density matrix in the vector form
# Output:
- the population in an array
"""
function get_population(rho_vec::Array{ComplexF64,1})
	N = size(rho_vec)[1]
	N = Int64(sqrt(N))
    P = zeros(Float64,N)
	# compute the trace for R
    for j = 1:N
        P[j] = P[j] + real(rho_vec[1+(j-1)*(N+1)])
    end
    return P
end

function get_population(rho_vec::Array{Float64,1})
	N = size(rho_vec)[1]
	N = Int64(sqrt(N))
    P = zeros(Float64,N)
	# compute the trace for R
    for j = 1:N
        P[j] = P[j] + rho_vec[1+(j-1)*(N+1)]
    end
    return P
end

function get_population_UQ_vector(rho_vec)
	N = size(rho_vec)[1]
	N = Int64(sqrt(N))
    #P = zeros(Float64,N)
	# compute the trace for R
	P = rho_vec[1]
    for j = 2:N
        P = [P;rho_vec[1+(j-1)*(N+1)]]
    end
    return P
end
"""
	get_population(RhoRealHistory::Array{Float64,2})

# Argument:
- RhoRealHistory: time history for the real part of the density matrix (saved as a flattened vector), j-th column corresponding to the j-th time point

rho_history = rhou_history - i rhov_history.

Since the density matrix should be Hermitian, it is good enough to just have the
real part.

# Output:
- The population history P which is a 2D array, P[i,j] is corresponding
to time point i and state j
"""
function get_population(rho_u_history::Array{Float64,2})
	# rho_u_history is RhoRealHistory
	N,nt = size(rho_u_history)
	N = Int64(sqrt(N))
    return transpose(rho_u_history[1:N+1:end,:])
end

function get_population(rho_u_history::Array{Array{Float64}})
	N,nt = size(rho_u_history)
	N = Int64(sqrt(N))
    return transpose(rho_u_history[1:N+1:end,:])
end

function get_population_UQ_matrix(rho_u_history)
	N,nt = size(rho_u_history)
	N = Int64(sqrt(N))
    return transpose(rho_u_history[1:N+1:end,:])
end

"""
	get_population(RhoVecHistory::Array{ComplexF64,2})

# Argument:
- RhoVecHistory, time history for the vectorized density matrix,
i-th column corresponding to the i-th time point
# Output:
- The population history P which is a 2D array, P[i,j] is corresponding
to time point i and state j
"""
function get_population(rho_vec_history::Array{ComplexF64,2})
	N,nt = size(rho_vec_history)
	N = Int64(sqrt(N))
	P = zeros(Float64,nt,N)
	# compute the trace
    return transpose(real(rho_vec_history[1:N+1:end,:]))
end

function get_population(rho_vec_history::Array{Array{ComplexF64}})
	N,nt = size(rho_vec_history)
	N = Int64(sqrt(N))
	P = zeros(Float64,nt,N)
	# compute the trace for R
	return transpose(real(rho_vec_history[1:N+1:end,:]))
end
"""
	convert_state_to_density(u::Array{Float64,1},v::Array{Float64,1})

# Argument:
- state vector: u-iv
# Output:
- vectorized density matrix: rho_u-i rho_v
"""
function convert_state_to_density(u::Array{Float64,1},v::Array{Float64,1})
	return (u*transpose(u)+v*transpose(v))[:],(v*transpose(u)-u*transpose(v))[:]
end
"""
	hello_world():
A hello world function
"""
function hello_world()
	println("Hello Quantum World!")
end
