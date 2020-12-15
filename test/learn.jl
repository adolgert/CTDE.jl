using Distributions
using Random
using Plots

rng = MersenneTwister(98274302)
wei = Weibull(3, 2)

x = zeros(100000)
rand!(rng, wei, x)
histogram(x)

trunc_wei = truncated(wei, 1.5, Inf)
y = zeros(100000)
rand!(rng, trunc_wei, y)
histogram(y)
