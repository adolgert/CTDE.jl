
"""
Weibull distribution. lambda is the scale parameter and k the
shape parameter.
F(T)=1-exp(-((T-Te)/lambda)^k)
"""
type TransitionWeibull <: TransitionDistribution
    parameters::Array{Float64,1}
end


function TransitionWeibull(lambda, k)
    TransitionWeibull([lambda, k, 0])
end


function TransitionWeibull(lambda, k, enabling_time)
    TransitionWeibull([lambda, k, enabling_time])
end


Parameters(tw::TransitionWeibull)=tw.parameters
function Parameters!(tw::TransitionWeibull, λ, k, Te)
    parameters[1]=λ
    parameters[2]=k
    parameters[3]=Te
end


function EnablingTime!(tw::TransitionWeibull, t::Float64)
    tw.parameters[3]=t
end


function Sample(dist::TransitionWeibull, now::Float64, rng::MersenneTwister)
    (λ, k, tₑ)=dist.parameters
    d=now-tₑ
    value=0
    U=Base.rand(rng)
    if d>0
        value=λ*(-log(1-U)+(d/λ)^k)^(1/k)-d
    else
        value=-d+λ*(-log(1-U))^(1/k)
    end
    now+value
end


function MeasuredSample(distribution::TransitionWeibull, now::Float64, rng)
    (λ, k, tₑ)=distribution.parameters
    d=now-tₑ
    value=0
    mlogU=randexp(rng)
    if d>0
        value=λ*(mlogU+(d/λ)^k)^(1/k)-d
    else
        value=-d+λ*(mlogU)^(1/k)
    end
    (now+value, mlogU)
end

function HazardIntegral(dist::TransitionWeibull, last, now)
    (λ, k, tₑ)=dist.parameters
    if now-tₑ>eps(Float64)
        return ((now-tₑ)/λ)^k - ((last-tₑ)/λ)^k
    else
        return 0::Float64
    end
end


function ConsumeSample(dist::TransitionWeibull, xa, start, finish)
    if xa<0
        xa=0
    end
    xa+HazardIntegral(dist, start, finish)
end


function CumulativeDistribution(dist::TransitionWeibull, when, now)
    1-exp(-hazard_integral(dist, now, when))
end


function ImplicitHazardIntegral(dist::TransitionWeibull,
        cumulative_hazard, when)
    (λ, k, tₑ)=dist.parameters
    if when-tₑ>eps(Float64)
        return tₑ + λ*(cumulative_hazard + ((when-tₑ)/λ)^k)^(1.0/k)
    else
        @debug("Weibull.implicit $cumulative_hazard")
        return tₑ + λ*(cumulative_hazard)^(1.0/k)
    end
end


function Putative(dist::TransitionWeibull, when,
        interval, consumed_interval)
    ImplicitHazardIntegral(dist, interval-consumed_interval, when)
end


function test(dist::TransitionWeibull)
    rng=MersenneTwister()
    (λ, k, tₑ)=dist.parameters
    ed=EmpiricalDistribution()
    for i in 1:10000
        push!(ed, Sample(dist, 0.0, rng))
    end
    expected_mean=λ*gamma(1+1/k)
    actual_mean=mean(ed)
    @debug("mean expected ", expected_mean, " actual ", actual_mean,
        " diff ", abs(expected_mean-actual_mean))

    expected_variance=λ^2*(gamma(1+2/k)-gamma(1+1/k)^2)
    obs_var=variance(ed)
    @debug("variance expected ", expected_variance, " actual ", obs_var,
        " diff ", abs(expected_variance-obs_var))

    min_value=min(ed)
    mink=min_value^k
    total=0.0
    for i in 1:length(ed)
        total+=ed.samples[i]^k
    end
    λ_estimator=(total/length(ed)-mink)^(1/k)
    @debug("λ expected ", λ, " actual ", λ_estimator,
        " diff ", abs(λ-λ_estimator))

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
    @debug("k expected ", k, " actual ", k_est,
        " diff ", abs(k-k_est))
end

