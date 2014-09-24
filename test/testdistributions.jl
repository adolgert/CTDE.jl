include("semimarkov.jl")

using SemiMarkov

function simple_test()
    exp_lambda=1.0
    exp_shift=0.0
    te=TransitionExponential(exp_lambda, exp_shift)
    test(TransitionExponential)
    wei_lambda=0.9
    wei_k=1.2
    wei_shift=0.0
    tw=TransitionWeibull(wei_lambda, wei_k, wei_shift)
    test(tw)
end

function shift_dist(dist, wrapped, tol)
    println("shift_dist ", typeof(dist))
    cnt=100
    for shift in [0.0, 0.2, 0.4]
        for i in 1:cnt
            U=i/cnt
            sample_time=wrapped.enabling_time+shift
            cdf_by_shift=quantile(wrapped, sample_time, U)
            cdf_by_hazard=implicit_hazard_integral(dist, -log(1-U), sample_time)
            if abs(cdf_by_shift-cdf_by_hazard)>tol
                println("U ", U, " Δ ", shift, " shift ", cdf_by_shift,
                    " λ ", cdf_by_hazard)
            else
                print(".")
            end
        end
    end
end


function shift_cdf(dist, wrapped, tol)
    println("shift_cdf ", typeof(dist))
    cnt=10
    for shift in [0.0, 0.2, 0.4]
        for i in 1:cnt
            sample_time=wrapped.enabling_time+shift
            T=sample_time+2.0*i/cnt
            cdf_by_shift=cdf(wrapped, T, sample_time)
            cdf_by_hazard=cdf(dist, T, sample_time)
            if abs(cdf_by_shift-cdf_by_hazard)>tol
                println("U ", U, " Δ ", shift, " shift ", cdf_by_shift,
                    " λ ", cdf_by_hazard)
            else
                print(".")
            end
        end
    end
    println()
end


function hazards_test()
    exp_lambda=1.0
    exp_shift=4.7
    te=TransitionExponential(exp_lambda, exp_shift)
    test(TransitionExponential)
    wei_lambda=0.9
    wei_k=1.2
    wei_shift=3.2
    tw=TransitionWeibull(wei_lambda, wei_k, wei_shift)
    test(tw)

    we=WrappedDistribution(Distributions.Exponential(1.0/exp_lambda), exp_shift)
    # shape k -> my k
    # scale s -> my lambda
    ww=WrappedDistribution(Distributions.Weibull(wei_k, wei_lambda), wei_shift)

    tol=1e-5
    println("exponential")
    shift_dist(te, we, tol)
    shift_cdf(te, we, tol)
    println("weibull")
    shift_dist(tw, ww, tol)
    shift_cdf(tw, ww, tol)
end

hazards_test()
