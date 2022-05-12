include("../../src/GLOQ.jl")
using LinearAlgebra
using Zygote
using Turing, Distributions, DifferentialEquations
# Import MCMCChain, Plots, and StatsPlots for visualizations and diagnostics.
using MCMCChains, Plots
using StatsPlots
using CSV,DataFrames
# Set a seed for reproducibility.
using Random
using DelimitedFiles
Random.seed!(16);
using Plots
using Dates
using StatsBase
using Distributions
using Printf
using KernelDensity

quad_strategy = 2

include("ReadDataAndDefineModel-BayesianRamsey-01-12-pm.jl")
# Read data
my_chain = read("chain-pm-2022-04-20/chain_500.jls",Chains)
my_chain_data = DataFrame(my_chain)
my_chain_size = 500
# Extract chains for essential values
freq_01_chain = my_chain_data._freq_01
freq_12_minus_chain = my_chain_data._freq_12_minus
freq_12_plus_chain = my_chain_data._freq_12_plus

freq_12_chain = [freq_12_minus_chain;freq_12_plus_chain]

freq_01_mean = mean(freq_01_chain)
freq_12_minus_mean = mean(freq_12_minus_chain)
freq_12_plus_mean = mean(freq_12_plus_chain)

# compute densities
@time kde_01    = kde(freq_01_chain,bandwidth=0.6e-7)
@time kde_plus  = kde(freq_12_plus_chain,bandwidth=1.5e-7)
@time kde_minus = kde(freq_12_minus_chain,bandwidth=1.5e-7)
@time kde_12 = kde(freq_12_chain,bandwidth=1e-7)

fig_01 = plot(kde_01,linewidth=2.5,
                size=(1200,800),
                xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
                legend=:false,
                title="Transition Frequency 0-1")
xticks!(fig_01,
        [freq_01_mean-6e-7;freq_01_mean;freq_01_mean+6e-7],
        ["mean - 0.6 KHz";string(@sprintf("%.6f",freq_01_mean)," GHz");"mean + 0.6 KHz"]);


fig_minus = plot(kde_minus,linewidth=2.5,
                    size=(1200,800),
                    xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
                    legend=:false,
                    title="Transition Frequency 1-2 minus parity flip")
xticks!(fig_minus,
        [freq_12_minus_mean-1.5e-6;freq_12_minus_mean;freq_12_minus_mean+1.5e-6],
        ["mean - 1.5 KHz";string(@sprintf("%.6f",freq_12_minus_mean)," GHz");"mean + 1.5 KHz"]);


fig_plus = plot(kde_plus,linewidth=2.5,
                   size=(1200,800),
                   xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
                   legend=:false,
                   title="Transition Frequency 1-2 plus parity flip")
xticks!(fig_plus,
       [freq_12_plus_mean-1.5e-6;freq_12_plus_mean;freq_12_plus_mean+1.5e-6],
       ["mean - 1.5 KHz";string(@sprintf("%.6f",freq_12_plus_mean)," GHz");"mean + 1.5 KHz"]);

display(fig_01)
display(fig_minus)
display(fig_plus)

freq_12_mean = mean(freq_12_chain)



fig_12 = plot(kde_12,line=2.5,
              size=(1200,800),
              title="Transition Frequency 1-2",
              xtickfontsize=18,ytickfontsize=18,titlefontsize=24,
              right_margin = 20Plots.mm,
              legend=:false)
xticks!(fig_12,
               [freq_12_mean-1.23e-4;freq_12_mean ;freq_12_mean+1.23e-4],
               ["mean - 1.23 KHz";string(@sprintf("%.6f",freq_12_mean)," GHz");"mean + 1.23 KHz"]);
display(fig_12)


# use spline and kde to construct a quadrature rule
#=
The following subroutine is written by Dr. Steven G. Johnson, see the notes
"Accurate solar-power integration: Solar-weighted Gaussian quadrature," arXiv:1912.06870 (2019),
Steven G. Johnson for more details.

The original code can be found in his Jupyter notebook: 
https://nbviewer.org/urls/math.mit.edu/~stevenj/Solar-Quadrature.ipynb

The authors are very grateful that Dr. Johnson shared his notebook with the community.
=#
using QuadGK, Dierckx
function gaussquad_interpolant(N::Integer,
                               X::AbstractVector{<:Real}, W::AbstractVector{<:Real},
                               a::Real=minimum(X), b::Real=maximum(X);
                               rtol::Real=sqrt(eps(typeof(float(b-a)))),
                               interpdegree::Integer=3)
    # some quick argument sanity checks:
    all(w ≥ 0 for w in W) || throw(ArgumentError("weights must be nonnegative"))
    length(X) == length(W) || throw(DimensionMismatch("points and weights must have same length"))
    N > 1 || throw(ArgumentError("a positive order N is required"))
    a < b || throw(ArgumentError("endpoints b > a are required"))
    rtol > 0 || throw(ArgumentError("a positive rtol is required"))
    interpdegree > 0 || throw(ArgumentError("a positive interpdegree is required"))
    
    # construct a cubic-spline interpolation from the data.
    # … we fit √W and then square it below to ensure non-negativity.
    Winterp_sqrt = Spline1D(X, sqrt.(W), k=interpdegree)
    
    # break integration interval at knots
    xab = sort!(collect(Iterators.filter(x -> a < x < b, Dierckx.get_knots(Winterp_sqrt))))
    push!(pushfirst!(xab, a), b) # add endpoints
    
    # quadrature routine for W-weighted integrands, that breaks integrand up into each
    # interpolation interval at the knots.
    interpquad(f, _a, _b; kws...) = quadgk(f, xab...; kws...)
    
    #println(rtol)
    return gauss(x -> Winterp_sqrt(x)^2, N, a, b; quad=interpquad, rtol=rtol)
    #return gauss(x -> Winterp_sqrt(x)^2, N, a, b; quad=interpquad, rtol=rtol)
end

N_quad_01 = 8
@time nodes_01,weights_01 = gaussquad_interpolant(N_quad_01,kde_01.x,kde_01.density,
                rtol=1e-10) 
nodes_01 .-= freq_01_mean

@time nodes_01_2,weights_01_2 = gauss(z->pdf(kde_01,z),N_quad_01,minimum(kde_01.x), maximum(kde_01.x));
nodes_01_2 .-= freq_01_mean

N_quad_12 = 8
@time nodes_plus,weights_plus = gaussquad_interpolant(N_quad_12,kde_plus.x,kde_plus.density,
                                rtol=1e-10) 
@time nodes_minus,weights_minus = gaussquad_interpolant(N_quad_12,kde_minus.x,kde_minus.density,
                                rtol=1e-10) 
nodes_12 = [nodes_minus;nodes_plus]
weights_12 = 0.5.*[weights_minus;weights_plus]

nodes_12 .-= freq_12_mean

# plot quadrature points
fig_nw_01 = plot(nodes_01,weights_01,
              size=(1200,800),
              linewidth=2.5,marker=(7,:circle),
              xlabel="Nodes",ylabel="Weights",
              title="Quadrature for 0-1",
              xtickfontsize=18,ytickfontsize=18,
              xlabelfontsize=18,ylabelfontsize=18,titlefontsize=24,
              left_margin = 5Plots.mm,bottom_margin= 5Plots.mm,
              legend=:false)
display(fig_nw_01)

fig_nw_12 = plot(nodes_12,weights_12,
              size=(1200,800),
              linewidth=2.5,marker=(7,:circle),
              xlabel="Nodes",ylabel="Weights",
              title="Quadrature for 1-2",
              xtickfontsize=18,ytickfontsize=18,
              xlabelfontsize=18,ylabelfontsize=18,titlefontsize=24,
              left_margin = 5Plots.mm,bottom_margin= 5Plots.mm,
              legend=:false)
display(fig_nw_12)


N_quad = N_quad_01*N_quad_12*2
nodes = zeros(N_quad,2)
weights = zeros(N_quad)
global ind = 1
for i = 1:N_quad_01
        for j = 1:N_quad_12*2
                global ind 
                nodes[ind,1] = nodes_01[i]
                nodes[ind,2] = nodes_12[j]
                weights[ind] = weights_01[i]*weights_12[j]
                ind += 1
        end
end

display(plot!(fig_nw_01,nodes_01_2,weights_01_2,marker=:square,line=:dash,legend=true))
# save the files
nodes_01_file   = string("nodes-weights/nodes01_",N_quad_01,".txt")
weights_01_file = string("nodes-weights/weights01_",N_quad_01,".txt")

nodes_12_file   = string("nodes-weights/nodes12_",N_quad_12*2,".txt")
weights_12_file = string("nodes-weights/weights12_",N_quad_12*2,".txt")

nodes_file   = string("nodes-weights/nodes-",N_quad_01,"_",N_quad_12*2,".txt")
weights_file = string("nodes-weights/weights-",N_quad_01,"_",N_quad_12*2,".txt")

writedlm(nodes_01_file,nodes_01);
writedlm(weights_01_file,weights_01);

writedlm(nodes_12_file,nodes_12);
writedlm(weights_12_file,weights_12);

writedlm(nodes_file,nodes);
writedlm(weights_file,weights);





