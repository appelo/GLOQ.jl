using Random
using Plots
using DifferentialEquations
pyplot()
#using GLOQ
#x = rand(Float64,8)
#X = [0 x[1]+x[2]*im x[3]+x[4]*im;
#     0 x[5]         x[6]+x[7]*im;
#     0 0            x[8]]
#H = 0.5*(X+adjoint(X))
H = [ 0 im 1;
     -im 0 0;
     1   0 0]
L1 = [0.0 0.25 0.0;
      0.0 0.0  0.5;
      0.0 0.0  0.0]
L2 = [0.0 0.0 0.0;
      0.0 0.5 0.0;
      0.0 0.0 1.0]
L_list = (L1,L2)
LH,LD = GLOQ.make_lindblad_operator(H,L_list)
L = LH+LD
println("Lindblad operator assembled")

tspan_length = 21
t_span = range(0.0,pi,length=tspan_length)
t_span2 = zeros(Float64,tspan_length)
dt = pi/(tspan_length-1);
for i = 1:tspan_length
      t_span2[i] = (i-1)*dt
end
u1 = [1.0;0;0]
rho1_history = GLOQ.exponential_solver(u1,L,t_span;initial_type = "states")
u2 = [0;1.0;0]
rho2 = (u2*u2')[:]
rho2_history = GLOQ.exponential_solver(rho2,L,t_span2;initial_type = "density")
u3 = [1.0/sqrt(2.0);-im/sqrt(2.0);0]
rho3_history = GLOQ.exponential_solver(u3,L,t_span;initial_type = "states")
u4 = [0.0;1.0/sqrt(2.0);-im/sqrt(2.0)]
rho4_history = GLOQ.exponential_solver(u4,L,t_span;initial_type = "states")
u5 = [1.0/sqrt(2.0);0.0;-im/sqrt(2.0)]
rho5_history = GLOQ.exponential_solver(u5,L,t_span;initial_type = "states")


population1_history = GLOQ.get_population(rho1_history)
population2_history = GLOQ.get_population(rho2_history)
population3_history = GLOQ.get_population(rho3_history)
population4_history = GLOQ.get_population(rho4_history)
population5_history = GLOQ.get_population(rho5_history)

fig = plot(t_span,population1_history,lab=:false)
plot!(fig,t_span,population2_history,lab=:false)
plot!(fig,t_span,population3_history,lab=:false)
plot!(fig,t_span,population4_history,lab=:false)
plot!(fig,t_span,population5_history,lab=:false)
display(fig)

# check the interface with DifferentialEquations
rho1 = (u1*u1')[:]
problem1 = GLOQ.LindbladODEProblem(rho1,L,t_span[end];initial_type="density")
rho1_history_diffeq = solve(problem1,saveat=t_span)
difference = maximum(abs.(rho1_history-rho1_history_diffeq))
@test(difference<5e-7)
