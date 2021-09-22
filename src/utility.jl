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

"""
	get_population(rho_vec_history::Array{ComplexF64,2})

# Argument:
- rho_vec_history, time history for the vectorized density matrix,
i-th column corresponding to the i-th time point
# Output:
- The population history P which is a 2D array, P[i,j] is corresponding
to time point i and state j
"""
function get_population(rho_vec_history::Array{ComplexF64,2})
	N,nt = size(rho_vec_history)
	N = Int64(sqrt(N))
	P = zeros(Float64,nt,N)
	# compute the trace for R
    for i = 1:nt
        for j = 1:N
            P[i,j] = P[i,j] + real(rho_vec_history[1+(j-1)*(N+1),i])
        end
    end
    return P
end
"""
	hello_world():
A hello world function
"""
function hello_world()
	println("Hello Quantum World!")
end
