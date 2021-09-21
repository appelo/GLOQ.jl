"""
  make_lindblad_operator(H,L_list,N::Int64=0)

  Argument:
  - the Hamiltonain H and a list of Lindblad terms L_k's
  - N is number of states

 Output:

 the operator for the vectorized system

       -LH: corresponding to the Hamiltonian part
       -LD: corresponding to the Dissipation part
"""
function make_lindblad_operator(H,L_list,N::Int64=0)
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
