The work flow for solving a quantum characterization problem consists of the following general steps:
1. Setup
2. Optimize

## 1. Setup
The setup phase includes specifying
- Load the experimental/synthetic data
- Specify the known parameters
- Number of the states
- Set up initial guess and the bounds for target parameters
- Create the loss (objective) function by fitting the data set with forward solves
- Optional: compute/approximate the Jacobian of the loss function. By default, auto-differentiation (Zygote.jl) or finite-difference approxiamtion will be used.

## 2. Optimization
- Feed the loss function and its Jacobian to the interface of optimization package.

