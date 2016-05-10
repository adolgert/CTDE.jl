
"""
A Gamma distribution.
α - shape parameter
β - inverse scale parameter, also called rate parameter

pdf=(1/(Γ(α)θ^α))x^(α-1) e^(-x/θ)
"""
function TransitionGamma(α::Float64, θ::Float64)
    # The supplied version uses θ=1/β.
    relative=Distributions.Gamma(α, θ)
    WrappedDistribution(relative)
end



"""
Given a distribution and an EmpiricalDistribution with
samples from it, compute a set of statistics on those samples,
returning a dictionary where the keys are the names of
the statistics and the entries are an array of
[expected, estimated].
"""
function TestGamma(d, ed::EmpiricalDistribution)
    statistic=Dict{AbstractString, Array}()
    (α, θ)=Parameters(d)
    expected_mean=αθ
    actual_mean=mean(ed)
    statistic["mean"]=[expected_mean, actual_mean]
    statistic["variance"]=[αθ^2, variance(ed)]

    min_value=min(ed)
    mink=min_value^k
    total=0.0
    for i in 1:length(ed)
        total+=ed.samples[i]^k
    end
    statistic["scale"]=[λ,
        (total/length(ed)-mink)^(1/k)]


    numerator=0.0
    denominator=0.0
    logsum=0.0
    for s in ed.samples
        numerator+=s^k*log(s) - mink*log(min_value)
        denominator+=s^k - mink
        logsum+=log(s)
    end
    k_est_inv=numerator/denominator - logsum/length(ed)
    k_est=1.0/k_est_inv

    statistic["exponent"]=[k, k_est]
    statistic
end


