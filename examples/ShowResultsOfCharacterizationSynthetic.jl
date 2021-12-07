########################
# Initial parameters
########################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 2.0*pi.*p_initial[1:3],omr_echo_01,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 0,
				 TC_01,t_echo_01,N_states)
initial_echo_01 = GLOQ.get_population(rho_Echo01_u)


# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 2.0*pi.*p_initial[1:3],omr_echo_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_echo_12,N_states)
initial_echo_12 = GLOQ.get_population(rho_Echo12_u)

# Echo 2-3
rho_Echo23_u,rho_Echo23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
				 2.0*pi.*p_initial[1:3],omr_echo_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_echo_23,N_states)
initial_echo_23 = GLOQ.get_population(rho_Echo23_u)

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 2.0*pi.*p_optim[1:3],omr_ramsey_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_ramsey_01,N_states)
initial_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)

# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 2.0*pi.*p_initial[1:3],omr_ramsey_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_ramsey_12,N_states)
initial_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)

# Ramsey 2-3
rho_Ramsey23_u,rho_Ramsey23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
				 2.0*pi.*p_initial[1:3],omr_ramsey_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_ramsey_23,N_states)
initial_ramsey_23 = GLOQ.get_population(rho_Ramsey23_u)


# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 2.0*pi.*p_initial[1:3],omr_t1_01,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 0,
				 TC_01,t_t1_01,N_states)
initial_t1_01 = GLOQ.get_population(rho_T101_u)

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 2.0*pi.*p_initial[1:3],omr_t1_12,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 1,
				 TC_12,t_t1_12,N_states)
initial_t1_12 = GLOQ.get_population(rho_T112_u)

# T1 2-3
rho_T123_u,rho_T123_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
				 2.0*pi.*p_initial[1:3],omr_t1_23,
				 2.0*pi.*p_initial[4:6],
				 p_initial[7:9],p_initial[10:12],
				 2,
				 TC_23,t_t1_23,N_states)
initial_t1_23 = GLOQ.get_population(rho_T123_u)

# Echo 0-1
initial_e01_error = plot(t_ramsey_01./1000.0,initial_echo_01-data_echo_01,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of Echo 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e01_error)

# Echo 1-2
initial_e12_error = plot(t_echo_12./1000.0,initial_echo_12-data_echo_12,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of Echo 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e12_error)

# Echo 2-3
initial_e23_error = plot(t_echo_23./1000.0,initial_echo_23-data_echo_23,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of Echo 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_e23_error)

# Ramsey 0-1
initial_r01_error = plot(t_ramsey_01./1000.0,initial_ramsey_01-data_ramsey_01,
		 size=(1200,600),legend=:outerright,
		 title = "Initial error of Ramsey 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r01_error)

# Ramsey 1-2
initial_r12_error = plot(t_ramsey_12./1000.0,initial_ramsey_12-data_ramsey_12,
				   size=(1200,600),legend=:outerright,
		   		   title="Initial error of Ramsey 1-2",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r12_error)

# Ramsey 2-3
initial_r23_error = plot(t_ramsey_23./1000.0,initial_ramsey_23-data_ramsey_23,
				   size=(1200,600),legend=:outerright,
		   		   title="Initial error of Ramsey 2-3",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_r23_error)

# T1 0-1
initial_t01_error = plot(t_t1_01./1000.0,initial_t1_01-data_t1_01,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of T1 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t01_error)

# T1 1-2
initial_t12_error = plot(t_t1_12./1000.0,initial_t1_12-data_t1_12,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of T1 1-2",
		 		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t12_error)

# T1 2-3
initial_t23_error = plot(t_t1_12./1000.0,initial_t1_23-data_t1_23,
		 size=(1200,600),legend=:outerright,
		 title = "Intial error of T1 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(initial_t23_error)

########################
# Optimized parameters
########################
# Echo 0-1
rho_Echo01_u,rho_Echo01_v = GLOQ.EchoParityForwardSolve(u0_echo_01,v0_echo_01,
				 2.0*pi.*p_optim[1:3],omr_echo_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_echo_01,N_states)
optim_echo_01 = GLOQ.get_population(rho_Echo01_u)


# Echo 1-2
rho_Echo12_u,rho_Echo12_v = GLOQ.EchoParityForwardSolve(u0_echo_12,v0_echo_12,
				 2.0*pi.*p_optim[1:3],omr_echo_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_echo_12,N_states)
optim_echo_12 = GLOQ.get_population(rho_Echo12_u)

# Echo 2-3
rho_Echo23_u,rho_Echo23_v = GLOQ.EchoParityForwardSolve(u0_echo_23,v0_echo_23,
				 2.0*pi.*p_optim[1:3],omr_echo_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_echo_23,N_states)
optim_echo_23 = GLOQ.get_population(rho_Echo23_u)

# Ramsey 0-1
rho_Ramsey01_u,rho_Ramsey01_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_01,v0_ramsey_01,
				 2.0*pi.*p_optim[1:3],omr_ramsey_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_ramsey_01,N_states)
optim_ramsey_01 = GLOQ.get_population(rho_Ramsey01_u)

# Ramsey 1-2
rho_Ramsey12_u,rho_Ramsey12_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_12,v0_ramsey_12,
				 2.0*pi.*p_optim[1:3],omr_ramsey_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_ramsey_12,N_states)
optim_ramsey_12 = GLOQ.get_population(rho_Ramsey12_u)

# Ramsey 2-3
rho_Ramsey23_u,rho_Ramsey23_v = GLOQ.RamseyParityForwardSolve(u0_ramsey_23,v0_ramsey_23,
				 2.0*pi.*p_optim[1:3],omr_ramsey_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_ramsey_23,N_states)
optim_ramsey_23 = GLOQ.get_population(rho_Ramsey23_u)


# T1 0-1
rho_T101_u,rho_T101_v = GLOQ.T1ParityForwardSolve(u0_t1_01,v0_t1_01,
				 2.0*pi.*p_optim[1:3],omr_t1_01,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 0,
				 TC_01,t_t1_01,N_states)
optim_t1_01 = GLOQ.get_population(rho_T101_u)

# T1 1-2
rho_T112_u,rho_T112_v = GLOQ.T1ParityForwardSolve(u0_t1_12,v0_t1_12,
				 2.0*pi.*p_optim[1:3],omr_t1_12,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 1,
				 TC_12,t_t1_12,N_states)
optim_t1_12 = GLOQ.get_population(rho_T112_u)

# T1 2-3
rho_T123_u,rho_T123_v = GLOQ.T1ParityForwardSolve(u0_t1_23,v0_t1_23,
				 2.0*pi.*p_optim[1:3],omr_t1_23,
				 2.0*pi.*p_optim[4:6],
				 p_optim[7:9],p_optim[10:12],
				 2,
				 TC_23,t_t1_23,N_states)
optim_t1_23 = GLOQ.get_population(rho_T123_u)



# Echo 0-1
fig_e01 = plot(t_echo_01./1000.0,data_echo_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"])
plot!(fig_e01,t_echo_01./1000.0,optim_echo_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim"],
          legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_e01)

fig_e01_error = plot(t_ramsey_01./1000.0,optim_echo_01-data_echo_01,
		 size=(1200,600),legend=:outerright,
		 title = "Error of Echo 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_e01_error)

# Echo 1-2
fig_e12 = plot(t_echo_12./1000.0,data_echo_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_e12,t_echo_12./1000.0,optim_echo_12,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_e12)

fig_e12_error = plot(t_echo_12./1000.0,optim_echo_12-data_echo_12,
		 size=(1200,600),legend=:outerright,
		 title = "Error of Echo 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_e12_error)

# Echo 2-3
fig_e23 = plot(t_echo_23./1000.0,data_echo_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Echo 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_e23,t_echo_23./1000.0,optim_echo_23,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_e23)

fig_e23_error = plot(t_echo_23./1000.0,optim_echo_23-data_echo_23,
		 size=(1200,600),legend=:outerright,
		 title = "Error of Echo 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_e23_error)

# Ramsey 0-1
fig_r01 = plot(t_ramsey_01./1000.0,data_ramsey_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_r01,t_ramsey_01./1000.0,optim_ramsey_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r01)

fig_r01_error = plot(t_ramsey_01./1000.0,optim_ramsey_01-data_ramsey_01,
		 size=(1200,600),legend=:outerright,
		 title = "Error of Ramsey 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_r01_error)

# Ramsey 1-2
fig_r12 = plot(t_ramsey_12./1000.0,data_ramsey_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_r12,t_ramsey_12./1000.0,optim_ramsey_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r12)

fig_r12_error = plot(t_ramsey_12./1000.0,optim_ramsey_12-data_ramsey_12,
				   size=(1200,600),legend=:outerright,
		   		   title="Error of Ramsey 1-2",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
                   legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_r12_error)

# Ramsey 2-3
fig_r23 = plot(t_ramsey_23./1000.0,data_ramsey_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="Ramsey 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_r23,t_ramsey_23./1000.0,optim_ramsey_23,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_r23)

fig_r23_error = plot(t_ramsey_23./1000.0,optim_ramsey_23-data_ramsey_23,
				   size=(1200,600),legend=:outerright,
		   		   title="Error of Ramsey 2-3",
		 	   	   label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
                   legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_r23_error)

# T1 0-1
fig_t01 = plot(t_t1_01./1000.0,data_t1_01,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 0-1",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_t01,t_t1_01./1000.0,optim_t1_01,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t01)

fig_t01_error = plot(t_t1_01./1000.0,optim_t1_01-data_t1_01,
		 size=(1200,600),legend=:outerright,
		 title = "Error of T1 0-1",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"])
display(fig_t01_error)

# T1 1-2
fig_t12 = plot(t_t1_12./1000.0,data_t1_12,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 1-2",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_t12,t_t1_12./1000.0,optim_t1_12,legend=:outerright,
		  label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t12)

fig_t12_error = plot(t_t1_12./1000.0,optim_t1_12-data_t1_12,
		 size=(1200,600),legend=:outerright,
		 title = "Error of T1 1-2",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_t12_error)

# T1 2-3
fig_t23 = plot(t_t1_23./1000.0,data_t1_23,line=(:dash),
		   size=(1200,600),legend=:outerright,
		   title="T1 2-3",
		   label=[L"$\rho_{00}$, data" L"$\rho_{11}$, data" L"$\rho_{22}$, data" L"$\rho_{33}$, data"],
           legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
plot!(fig_t23,t_t1_23./1000.0,optim_t1_23,legend=:outerright,label=[L"$\rho_{00}$, sim" L"$\rho_{11}$, sim" L"$\rho_{22}$, sim" L"$\rho_{33}$, sim"])
display(fig_t23)

fig_t23_error = plot(t_t1_23./1000.0,optim_t1_23-data_t1_23,
		 size=(1200,600),legend=:outerright,
		 title = "Error of T1 2-3",
		 label=[L"$\rho_{00}$" L"$\rho_{11}$" L"$\rho_{22}$" L"$\rho_{33}$"],
         legendfontsize=15,xtickfontsize=15,ytickfontsize=15,titlefontsize=18)
display(fig_t23_error)
