"""
    make_lindblad_operator(HK,HS,L_list,N::Int64=0)

# Argument:
- The real part and the imaginary part of the Hamiltonian: ``H = H_K-i H_S``
- A list of Lindblad terms L_k's
- N is number of states

# Output:
the operator for the vectorized system

- LK: real matrix corresponding to the real part of the Hamiltonian operator
- LS: real matrix corresponding to the imaginary part of the Hamiltonian operator
- LD: real matrix corresponding to the Lindblad operator
- Final system: ``\\rho_t = -i(L_K+i L_S)\\rho+L_D \\rho``
"""
function make_lindblad_operator(HK::Array{Float64,2},HS::Array{Float64,2},L_list,N::Int64=0)
    if(N==0)
        N = size(HK)[1]
    end
    # define an identity operator
    It = Array{Float64, 2}(I, N, N)
    # Hamiltonian part
    LK = (kron(It,HK) - kron(transpose(HK),It))
    LS = (kron(It,HS) - kron(transpose(HS),It))
    # Disspation part
    LD = zeros(Float64,N^2,N^2)
    for k = 1:length(L_list)
        LD += kron(L_list[k],L_list[k])-
              0.5*(kron(It,transpose(L_list[k])*L_list[k] )+kron(transpose(transpose(L_list[k])*L_list[k]),It))
    end
    return LK,LS,LD
end

"""
    make_lindblad_operator(H,L_list,N::Int64=0)

# Argument:
- the Hamiltonain H and a list of Lindblad terms L_k's
- N is number of states

# Output:
the operator for the vectorized system

- LH: complex matrix corresponding to the Hamiltonian operator
- LD: complex matrix corresponding to the Lindblad operator
"""
function make_lindblad_operator_complex(H,L_list,N::Int64=0)
    if(N==0)
        N = size(H)[1]
    end
    # define an identity operator
    It = Array{Float64, 2}(I, N, N)
    # Hamiltonian part
    LH = -1.0im*(kron(It,H) - kron(transpose(H),It))
    # Disspation part
    LD = zeros(Float64,N^2,N^2)
    for k = 1:length(L_list)
        LD += kron(L_list[k],L_list[k])-
              0.5*(kron(It,transpose(L_list[k])*L_list[k] )+kron(transpose(transpose(L_list[k])*L_list[k]),It))
    end
    return LH,LD
end
