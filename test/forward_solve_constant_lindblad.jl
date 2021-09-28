#using GLOQ
#x = rand(Float64,8)
#X = [0 x[1]+x[2]*im x[3]+x[4]*im;
#     0 x[5]         x[6]+x[7]*im;
#     0 0            x[8]]
#H = 0.5*(X+adjoint(X))
HK = [ 0.0 0.0 1.0;
       0.0 0.0 0.0;
       1.0 0.0 0.0]
HS = [ 0.0 1.0 0.0;
      -1.0 0.0 0.0;
       0.0 0.0 0.0]
H = HK + 1.0im.*HS
L1 = [0.0 0.25 0.0;
      0.0 0.0  0.5;
      0.0 0.0  0.0]
L2 = [0.0 0.0 0.0;
      0.0 0.5 0.0;
      0.0 0.0 1.0]
L_list = (L1,L2)
LH,LD = GLOQ.make_lindblad_operator_complex(H,L_list)
LK,LS,LD2 = GLOQ.make_lindblad_operator(HK,HS,L_list)
L = LH+LD
println("Test whether the real-valued and complex-valued operator match each other")
@test (LD==LD2)&&( LH==-1.0im*(LK+1.0im*LS) )

println("\n Lindblad operator assembled\n")

tspan_length = 21
t_span = range(0.0,pi,length=tspan_length)
t_span2 = zeros(Float64,tspan_length)
dt = pi/(tspan_length-1);
for i = 1:tspan_length
      t_span2[i] = (i-1)*dt
end
u1 = [1.0;0;0]
rho1_history = GLOQ.exponential_solver_complex(u1,L,t_span;initial_type = "states")
u2 = [0;1.0;0]
rho2 = (u2*u2')[:]
rho2_history = GLOQ.exponential_solver_complex(rho2,L,t_span2;initial_type = "density")
u3 = [1.0/sqrt(2.0);-im/sqrt(2.0);0]
rho3_history = GLOQ.exponential_solver_complex(u3,L,t_span;initial_type = "states")
u4 = [0.0;1.0/sqrt(2.0);-im/sqrt(2.0)]
rho4_history = GLOQ.exponential_solver_complex(u4,L,t_span;initial_type = "states")
u5 = [1.0/sqrt(2.0);0.0;-im/sqrt(2.0)]
rho5_history = GLOQ.exponential_solver_complex(u5,L,t_span;initial_type = "states")


population1_history = GLOQ.get_population(rho1_history)
population2_history = GLOQ.get_population(rho2_history)
population3_history = GLOQ.get_population(rho3_history)
population4_history = GLOQ.get_population(rho4_history)
population5_history = GLOQ.get_population(rho5_history)

#################################################
# test real solvers
rho1_history_u,rho1_history_v = GLOQ.exponential_solver(real(u1),-imag(u1),LK,LS,LD2,t_span;initial_type = "states")
rho2_u,rho2_v = GLOQ.convert_state_to_density(real(u2),-imag(u2))
rho2_history_u,rho2_history_v = GLOQ.exponential_solver(rho2_u,rho2_v,LK,LS,LD2,t_span;initial_type = "density")
rho3_history_u,rho3_history_v = GLOQ.exponential_solver(real(u3),-imag(u3),-imag(LH),real(LH),LD,t_span;initial_type = "states")
rho4_u,rho4_v = GLOQ.convert_state_to_density(real(u4),-imag(u4))
rho4_history_u,rho4_history_v = GLOQ.exponential_solver(rho4_u,rho4_v,LK,LS,LD2,t_span;initial_type = "density")
rho5_history_u,rho5_history_v = GLOQ.exponential_solver(real(u5),-imag(u5),LK,LS,LD2,t_span;initial_type = "states")

@test ( (norm(rho1_history_u-1.0im*rho1_history_v-rho1_history,Inf)<1e-14)&&
        (norm(rho2_history_u-1.0im*rho2_history_v-rho2_history,Inf)<1e-14)&&
        (norm(rho3_history_u-1.0im*rho3_history_v-rho3_history,Inf)<1e-14)&&
        (norm(rho4_history_u-1.0im*rho4_history_v-rho4_history,Inf)<1e-14)&&
        (norm(rho5_history_u-1.0im*rho5_history_v-rho5_history,Inf)<1e-14) )

population1_history2 = GLOQ.get_population(rho1_history_u)
population2_history2 = GLOQ.get_population(rho2_history_u)
population3_history2 = GLOQ.get_population(rho3_history_u)
population4_history2 = GLOQ.get_population(rho4_history_u)
population5_history2 = GLOQ.get_population(rho5_history_u)

@test ( (norm(population1_history-population1_history2,Inf)<1e-14)&&
        (norm(population2_history-population2_history2,Inf)<1e-14)&&
        (norm(population3_history-population3_history2,Inf)<1e-14)&&
        (norm(population4_history-population4_history2,Inf)<1e-14)&&
        (norm(population5_history-population5_history2,Inf)<1e-14) )


# check the interface with DifferentialEquations

rho3 = (u3*u3')[:]
problem3 = GLOQ.LindbladODEProblemComplex(rho3,L,t_span[end];initial_type="density")
rho3_history_diffeq = solve(problem3,saveat=t_span)
difference = maximum(abs.(rho3_history-rho3_history_diffeq))
@test(difference<5e-7)

N_density = length(u3)^2;
problem3_real = GLOQ.LindbladODEProblem(real(u3),-imag(u3),LK,LS,LD,t_span[end];initial_type="states")
rho3_history_diffeq_real = solve(problem3_real,saveat=t_span)

println(maximum(abs.(real(rho3_history)-rho3_history_diffeq_real[1:N_density,:])))
println(maximum(abs.(imag(rho3_history)+rho3_history_diffeq_real[N_density+1:end,:])))
@test( (maximum(abs.(real(rho3_history)-rho3_history_diffeq_real[1:N_density,:]))<5e-5)&&
       (maximum(abs.(imag(rho3_history)+rho3_history_diffeq_real[N_density+1:end,:]))<5e-5) )
