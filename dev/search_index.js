var documenterSearchIndex = {"docs":
[{"location":"methods/","page":"Methods","title":"Methods","text":"The following methods (functions) are exported and available by using GLOQ.","category":"page"},{"location":"methods/","page":"Methods","title":"Methods","text":"Modules = [GLOQ]\nOrder = [:function]","category":"page"},{"location":"methods/#GLOQ.LindbladODEProblem-Tuple{Any,Array{Complex{Float64},2},Float64}","page":"Methods","title":"GLOQ.LindbladODEProblem","text":"LindbladODEProblem(rho0,L::Array{ComplexF64,2},time_final::Float64;initial_type = \"density\")\n\nFunction provide interfaces to DifferentialEquations package.\n\nArgument:\n\nL: Lindblad operator\ntime: final time or time beining evaluated\nrho0: initial condition\ninitial_type: specify initial value is a density matrix/a state vector\n\nOutput:\n\nA problem object which we will feed to DifferentialEquations.jl\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.about-Tuple{}","page":"Methods","title":"GLOQ.about","text":"About GLOQ.jl\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.exponential_solver-Tuple{Any,Any,Array{Float64,N} where N}","page":"Methods","title":"GLOQ.exponential_solver","text":"exponential_solver(rho_vec0,L,t_span::Array{Float64};initial_type = \"density\"):\n\nArgument:\n\nL: the whole propagation operator\nt_span: where the function will be evaluated stored in\nrho0_vec: initial density/state\ninitial_type: decide we aer given an initial density matrix or a state vector\n\nOutput:\n\nsolutions at t_span\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.exponential_solver-Tuple{Any,Any,StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}}","page":"Methods","title":"GLOQ.exponential_solver","text":"exponential_solver(rho_vec0,L,t_span::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}};\n                   initial_type = \"density\")\n\nArgument:\n\nL: the whole propagation operator\nt_span: where the function will be evaluated stored in\nrho0_vec: initial density/state\ninitial_type: decide we aer given an initial density matrix or a state vector\n\nOutput:\n\nsolutions at t_span\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.get_population-Tuple{Array{Complex{Float64},1}}","page":"Methods","title":"GLOQ.get_population","text":"get_population(rho_vec::Array{ComplexF64,1})\n\nArgument:\n\nrho_vec, density matrix in the vector form\n\nOutput:\n\nthe population in an array\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.get_population-Tuple{Array{Complex{Float64},2}}","page":"Methods","title":"GLOQ.get_population","text":"get_population(rho_vec_history::Array{ComplexF64,2})\n\nArgument:\n\nrhovechistory, time history for the vectorized density matrix,\n\ni-th column corresponding to the i-th time point\n\nOutput:\n\nThe population history P which is a 2D array, P[i,j] is corresponding\n\nto time point i and state j\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.hello_world-Tuple{}","page":"Methods","title":"GLOQ.hello_world","text":"hello_world():\n\nA hello world function\n\n\n\n\n\n","category":"method"},{"location":"methods/#GLOQ.make_lindblad_operator","page":"Methods","title":"GLOQ.make_lindblad_operator","text":"make_lindblad_operator(H,L_list,N::Int64=0)\n\nArgument:\n\nthe Hamiltonain H and a list of Lindblad terms L_k's\nN is number of states\n\nOutput:\n\nthe operator for the vectorized system\n\nLH: corresponding to the Hamiltonian part\nLD: corresponding to the Dissipation part\n\n\n\n\n\n","category":"function"},{"location":"#GLOQ.jl","page":"Home","title":"GLOQ.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for GLOQ.jl","category":"page"},{"location":"function-index/","page":"Index","title":"Index","text":"Modules = [GLOQ]","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"The following types are exported and available by using GLOQ.","category":"page"},{"location":"types/","page":"Types","title":"Types","text":"Modules = [GLOQ]\nOrder = [:type]","category":"page"}]
}
