The `GLOQ.jl` package aims at solving quantum characterization problems with deterministic optimization or Bayesian inference methods. 


The general work flow for solving a quantum characterization problem consists of the following general steps:
1. Setup the characterization problem.
2. Perform optimization/Bayesian inference.

## Deterministic optimization
### 1. Setup
The setup phase includes 
- Loading the experimental/synthetic data.
- Specify the known parameters.
- Specify the number of the states of the device.


### 2. Optimization
- Set up an initial guess and bounds for the parameters to be characterized. 
- Create the objective function by fitting the data set with forward solves
- Optional: compute/approximate the Jacobian of the loss function. By default, auto-differentiation (Zygote.jl) or finite-difference approximation will be used.
- Feed the loss function and its Jacobian to the interface of optimization package.

## Bayesian inference with Turing.jl
### 1. Setup
The setup phase includes specifying
- Load the experimental/synthetic data
- Specify the known parameters
- Number of the states


### 2. Bayesian inference with [Turing.jl](https://turing.ml/stable/)
- Use `@model` interface of `Turing.jl` to define a `model` object fitting the data based on some random parameter values sampled from some priori and proposal distributions.
- Feed the `model` object to the interface of a Markov Chain Monte-Carlo (MCMC) sampler.
