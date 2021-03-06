# GLOQ.jl

GLOQ.jl (Göran Lindblad Open Quantum systems) is a Julia package solving characterization problems for open quantum systems based on Lindblad equations with
deterministic optimization method and Bayesian inference.

GLOQ.jl provides Julia native implementations of forward Lindblad solves with interfaces to the differential equation solver [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/), a unified Julia interface to various optimization packages [GalacticOptim.jl](https://galacticoptim.sciml.ai/stable/), and the Julia Bayesian inference package [Turing.jl](https://turing.ml/stable/) (partially supported for now). 

