# Lindblad equation
We consider the $N$ energy level Lindblad equation:
```math
\begin{equation}
\dot{\rho} = -i\left(H\rho - \rho H\right) + \sum_{j=1}^2 \left( {\cal L}_{j} \rho {\cal L}_{j}^\dagger -
\frac{1}{2}\left( {\cal L}_{j}^\dagger{\cal L}_{j}\rho + \rho{\cal L}_{j}^\dagger{\cal L}_{j} \right) \right).
\end{equation}
```
Here, $\rho(t)$ is the density matrix, $H(t)$ is a Hamiltonain operator, and operators ${\cal L}_{1}$ and ${\cal L}_{2}$ model the decay and dephasing mechanisim.

By default, the Hamiltonian is in the form of $H(t)=H_{\rm 0}+H_c(t)$, where
```math
H_{\rm 0} = \left[
\begin{array}{cccc}
0 &  &  &   & \\
  & \omega_{01} &  & & \\
  &   & \omega_{01} + \omega_{12}  &  &\\
  &   &  &  \ddots  & \\
  &   &  &          & \sum_{k=0}^{N-1}\omega_{N-1,N}
\end{array}
\right]
```
is the system Hamiltonian, and $H_c(t)$ is the control Hamiltonian.
The decay operator ${\cal L}_{1}$ is in the form 
```math
{\cal L}_{1} = \left[
\begin{array}{cccc}
0 & \sqrt{\gamma_{1,1}} &  &   &\\
  & 0 & \sqrt{\gamma_{1,2}}  &  &\\
  &   & \ddots & \ddots &  \\
  &   &   &   & 0  &\sqrt{\gamma_{1,N-1}} \\
  &   &   &   &    & 0
\end{array}\right],
```
and the dephasing operator ${\cal L}_{2}$ is in the form 
```math
{\cal L}_{2} = \left[
\begin{array}{cccc}
0 &   &   &  & \\
  & \sqrt{\gamma_{2,1}} &   &   & \\
  &   & \sqrt{\gamma_{2,2}}  &  & \\
  &   &   &  \ddots  & \\
  &   &   &          & \sqrt{\gamma_{2,N-1}}
\end{array}
\right].
```
The $T_1$ relaxation time for energy level $k$ is $T_{1,k}=1/\gamma_{1,k}$,
and the $T_2$ dephasing time for energy level $k$ is $T_{2,k}=1/\gamma_{2,k}$.

## Lindblad equation with charge noise
We also consider the Lindblad model with noise. The system Hamiltonian is
```math
H_{\rm 0} = \left[
\begin{array}{cccc}
0 &  &  &   & \\
  & \widetilde{\omega}_{01} &  & & \\
  &   & \widetilde{\omega}_{01} + \widetilde{\omega}_{12}  &  &\\
  &   &  &  \ddots  & \\
  &   &  &          & \sum_{k=0}^{N-1}\widetilde{\omega}_{N-1,N}
\end{array}
\right],
```
and 
$$\widetilde{\omega}_{k,k+1}=\omega_{k,k+1}-\pi d_{k,k+1}\cos\left(\pi\frac{C+pe}{e}\right).$$
Here, $p\in\{0,1\}$ is the parity, $d_{k,k+1}$ is the charge dispersion,
$C$ is the charge over the Josephson junction, $e$ is the charge of an electron, and $\frac{C}{e}\in[0,1]$.

For now, we neglect the charge noise due to the charge event (the effect of $C/e$) and only consider the charge noise due to the parity event (the effect of $p$). 
In the numerical simulation, we take the average of the forward solve of the deterministic Lindblad equation with $p=0$ and $p=1$. 