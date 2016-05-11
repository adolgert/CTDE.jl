using UNURAN
###################################################
# These functions help with any unuran distribution
###################################################

function MakeExponential()
	distr=unur_distr_exponential(C_NULL, 0)
	par=unur_auto_new(distr)
	gen=unur_init(par)
	unur_distr_free(distr)
	gen
end

g_Unit_Exponential=MakeExponential()

function TruncateDomain(gen::Ptr{UNUR_GEN}, left)
	distr=unur_get_distr(gen)
	unur_distr_cont_set_domain(distr, left, Inf)
	assert(unur_reinit(gen)==UNUR_SUCCESS)
end

function TruncateDiscrDomain(gen::Ptr{UNUR_GEN}, times, left)
    # Times are the times of the jumps in the pmf.
    # They are strictly increasing.
    distr=unur_get_distr(gen)
    ileft=searchsortedfirst(times, left)
    unur_distr_discr_set_domain(distr, ileft, Inf)
    assert(unur_reinit(gen)==UNUR_SUCCESS)
end

function ResetDomain(gen::Ptr{UNUR_GEN})
	distr=unur_get_distr(gen)
	unur_distr_cont_set_domain(distr, -Inf, Inf)
	assert(unur_reinit(gen)==UNUR_SUCCESS)
end

# Two different ways to sample and get the quantile.
# This one samples, then asks for the quantile.
function SampleThenQuantile(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR})
	v=unur_sample_cont(gen)
	cdf=unur_distr_cont_get_cdf(distr)
	S=1.0-ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}), v, distr)
	(v, S)
end

# This one samples, then asks for the quantile.
function SampleDiscrThenQuantile(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR})
    v=unur_sample_discr(gen)
    S=1-unur_distr_discr_eval_cdf(v, distr)
    (v, S)
end

# This one samples by inversion.
function SampleByInversion(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR},
		rng::Ptr{UNUR_URNG})
	U=unur_urng_sample(rng)
	v=unur_quantile(gen, U)
	(v, 1-U)
end

function ConditionalSurvival(distr::Ptr{UNUR_DISTR}, a::Int, b::Int)
	cdf=unur_distr_cont_get_cdf(distr)
	if a>eps(Float64)
    	last=ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}),
    			a, distr)
    else
    	last=0.0
    end
	current=ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}),
			b, distr)
	(1.0-last)/(1.0-current)
end


function ConditionalDiscrSurvival(distr::Ptr{UNUR_DISTR}, a, b)
    if a>eps(Float64)
        last=unur_distr_discr_eval_cdf(a, distr)
    else
        last=0.0
    end
    current=unur_distr_discr_eval_cdf(b, distr)
    (1.0-last)/(1.0-current)
end

# We are given the conditional survival for a distribution
# from some time (called shift) to a future time.
# 1. Find the survival of the distribution to that shifted time.
# 2. Use that and the conditonal survival to invert the distribution.
function ShiftedQuantile(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR},
		survival, shift)
	s_prime=unur_distr_discr_eval_cdf(shift, distr)
  	unur_quantile(gen, 1.0-survival*s_prime)
end

function randexp(rng::Ptr{UNUR_URNG})
	unur_sample_cont(g_Unit_Exponential::Ptr{UNUR_GEN})
end

function rand(rng::Ptr{UNUR_URNG})
    unur_sample_urng(rng)
end


####################################################
# Continuous Unur Distributions can derive from this
####################################################

abstract UnurDistribution <: TransitionDistribution


function EnablingTime!(gu::UnurDistribution, t::Float64)
    gu.te=t
end


function Sample(dist::UnurDistribution, now::Float64, rng)
    d=now-dist.te
    if abs(d)<eps(1.0)
        v=dist.te+unur_sample_cont(dist.gen)
    else
        TruncateDomain(dist.gen, d)
        v=dist.te+unur_sample_cont(dist.gen)
    end
    v
end


function MeasuredSample(gu::UnurDistribution, now::Float64, rng)
    d=now-gu.te
    if abs(d)<eps(1.0)
        (v, S)=SampleThenQuantile(gu.gen, gu.distr)
    else
        TruncateDomain(gu.gen, d)
        (v, S)=SampleThenQuantile(gu.gen, gu.distr)
    end
    (gu.te+v, S)
end


# xa is the conditonal survival
function ConsumeSample(gu::UnurDistribution, xa, start, finish)
    xa=(xa<0) ? 1 : xa
    xa*ConditionalSurvival(gu.distr, start-gu.te, finish-gu.te)
end


function Putative(gu::UnurDistribution, when, interval, consumed_interval)
    # If you want unur_quantile to work, you have to use
    # HINV, NINV, PINV, CSTD, or DGT. CSTD requires that
    # the generator have an inversion method.
    S=interval/consumed_interval
    shift=when-gu.te
    gu.te+ShiftedQuantile(gu.gen, gu.distr, S, shift)
end


##########################################
# Gamma Distribution
##########################################
export GammaUnur

# Gamma Distribution, implemented with UNURAND
# http://statmath.wu.ac.at/unuran/doc/unuran.html#gamma
type GammaUnur <: UnurDistribution
	parameters::Array{Float64,1}
	te::Float64
	gen::Ptr{UNUR_GEN}
	distr::Ptr{UNUR_DISTR}
end

function MakeGammaUnur(rng, alpha, beta)
    distr=unur_distr_gamma([alpha, beta], 2)
    # The auto method uses transformed density rejection,
    # which works great when alpha>1, which it isn't
    # for the problem at hand.
    par=unur_ninv_new(distr)
    unur_set_urng(par, rng)
    gen=unur_init(par)
    (gen, distr)
end

function GammaUnur(rng::Ptr{UNUR_URNG}, alpha, beta)
    gen, distr=MakeGammaUnur(rng, alpha, beta)
	gu=GammaUnur([alpha, beta], 0.0, gen, distr)
	gu
end


Parameters(gu::GammaUnur)=gu.parameters
function Parameters!(gu::GammaUnur, rng, alpha, beta)
	gu.parameters[1]=alpha
	gu.parameters[2]=beta
	if gu.gen!=nothing
		unur_free(gu.gen)
		unur_distr_free(gu.distr)
	end
    gu.gen, gu.distr=MakeGammaUnur(rng, alpha, beta)
end


##############################################
# Uniform distribution
##############################################


export UniformUnur

type UniformUnur <: UnurDistribution
    parameters::Array{Float64,1}
    te::Float64
    gen::Ptr{UNUR_GEN}
    distr::Ptr{UNUR_DISTR}
end


function UniformUnur(rng::Ptr{UNUR_URNG}, alpha, beta)
    distr=unur_distr_uniform([alpha, beta], 2)
    par=unur_auto_new(distr)
    unur_set_urng(par, rng)
    gen=unur_init(par)
    gu=UniformUnur([alpha, beta], 0.0, gen, distr)
    gu
end


Parameters(gu::UniformUnur)=gu.parameters
function Parameters!(gu::UniformUnur, rng, alpha, beta)
    gu.parameters[1]=alpha
    gu.parameters[2]=beta
    if gu.gen!=nothing
        unur_free(gu.gen)
        unur_distr_free(gu.distr)
    end
    gu.distr=unur_distr_uniform([alpha, beta], 2)
    out_params=Ref{Ptr{Cdouble}}(0)
    nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
    arr=pointer_to_array(out_params[], nparams)
    par=unur_auto_new(gu.distr)
    unur_set_urng(par, rng)
    gu.gen=unur_init(par)
end


#########################################
# Log-logistic distribution
#########################################

export LogLogisticUnur

function LogLogisticPdf(x, distr::Ptr{UNUR_DISTR})
    out_params=Ref{Ptr{Cdouble}}(0)
    nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
    arr=pointer_to_array(out_params[], nparams)
    a=arr[1]
    b=arr[2]
    v=(b/a)*(x/a)^(b-1)/( 1+(x/a)^b )^2
    convert(Cdouble, v)::Cdouble
end

function LogLogisticDPdf(x, distr::Ptr{UNUR_DISTR})
    out_params=Ref{Ptr{Cdouble}}(0)
    nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
    arr=pointer_to_array(out_params[], nparams)
    a=arr[1]
    b=arr[2]
    v=(a^(1-2*b)*x^(b-2) *((b-1)*a^b *((a/x)^(-b)+1)-2*b*x^b))/((a/x)^(-b)+1)^3
    convert(Cdouble, v)::Cdouble
end

function LogLogisticCdf(x, distr::Ptr{UNUR_DISTR})
    out_params=Ref{Ptr{Cdouble}}(0)
    nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
    arr=pointer_to_array(out_params[], nparams)
    a=arr[1]
    b=arr[2]
    v=1.0/( 1+(x/a)^b )
    convert(Cdouble, v)::Cdouble
end

# y=1/(1+(x/a)^b)
# 1+(x/a)^b = 1/y
# (x/a)^b = 1/y - 1
# (x/a)=((1/y)-1)^(1/b)
# x=a((1/y)-1)^(1/b)
function LogLogisticInvCdf(x, distr::Ptr{UNUR_DISTR})
    out_params=Ref{Ptr{Cdouble}}(0)
    nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
    arr=pointer_to_array(out_params[], nparams)
    a=arr[1]
    b=arr[2]
    v=a*((1/y)-1)^(1/b)
    convert(Cdouble, v)::Cdouble
end

function MakeLogLogistic(a, b)
    const pdf=cfunction(LogLogisticPdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}))
    const dpdf=cfunction(LogLogisticDPdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}))
    const cdf=cfunction(LogLogisticCdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}))
    const inv=cfunction(LogLogisticInvCdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}))
    distr=unur_distr_cont_new()
    unur_distr_cont_set_pdf( distr, pdf )
    unur_distr_cont_set_dpdf( distr, dpdf )
    unur_distr_cont_set_cdf( distr, cdf )
    unur_distr_cont_set_invcdf( distr, inv )
    unur_distr_cont_set_domain( distr, 0, UNUR_INFINITY)
    if b>1
        unur_distr_cont_set_mode(distr, a*( (b-1)/(b+1) )^(1/b))
    else
        unur_distr_cont_set_mode(distr, 0)
    end
    par=unur_tdr_new(distr)
    gen=unur_init(par)
    (distr, gen)
end


export LogLogisticUnur

# PDF is (b/a)(x/a)^(b-1)/( 1+(x/a)^b )^2
# CDF is 1/( 1+(x/a)^(-b) )
type LogLogisticUnur <: UnurDistribution
    parameters::Array{Float64,1}
    te::Float64
    gen::Ptr{UNUR_GEN}
    distr::Ptr{UNUR_DISTR}
end


function LogLogisticUnur(rng::Ptr{UNUR_URNG}, alpha, beta)
    distr, gen=MakeLogLogistic(alpha, beta)
    LogLogisticUnur([alpha, beta], 0.0, gen, distr)
end


Parameters(gu::LogLogisticUnur)=gu.parameters
function Parameters!(gu::LogLogisticUnur, rng, alpha, beta)
    gu.parameters[1]=alpha
    gu.parameters[2]=beta
    if gu.gen!=nothing
        unur_free(gu.gen)
        unur_distr_free(gu.distr)
    end
    distr, gen=MakeLogLogistic(alpha, beta)
    LogLogisticUnur([alpha, beta], 0.0, gen, distr)
end


########################
# Discrete distribution
########################

abstract UnurDiscrDistribution <: TransitionDistribution


function EnablingTime!(gu::UnurDiscrDistribution, t::Float64)
    gu.te=t
end


function Sample(dist::UnurDiscrDistribution, now::Float64, rng)
    d=now-dist.te
    if abs(d)<eps(1.0)
        v=dist.te+gu.times[unur_sample_discr(dist.gen)]
    else
        TruncateDiscrDomain(dist.gen, gu.times, d)
        v=dist.te+gu.times[unur_sample_cont(dist.gen)]
    end
    v
end


function MeasuredSample(gu::UnurDiscrDistribution, now::Float64, rng)
    d=now-gu.te
    if abs(d)<eps(1.0)
        (v, S)=SampleDiscrThenQuantile(gu.gen, gu.distr)
    else
        TruncateDiscrDomain(gu.gen, gu.times, d)
        (v, S)=SampleDiscrThenQuantile(gu.gen, gu.distr)
    end
    (gu.te+gu.times[v], S)
end


# xa is the conditonal survival
function ConsumeSample(gu::UnurDiscrDistribution, xa, start, finish)
    xa=(xa<0) ? 1 : xa
    xa*ConditionalDiscrSurvival(gu.distr, start-gu.te, finish-gu.te)
end


function Putative(gu::UnurDiscrDistribution, when, interval, consumed_interval)
    # If you want unur_quantile to work, you have to use
    # HINV, NINV, PINV, CSTD, or DGT. CSTD requires that
    # the generator have an inversion method.
    S=interval/consumed_interval
    shift=when-gu.te
    gu.te+ShiftedQuantile(gu.gen, gu.distr, S, shift)
end



# This distribution is made by specifying the probability
# mass function as an array.
function MakeDiscretePVUnur(rng, pmf::Array{Float64,1})
    distr=unur_distr_discr_new()
    unur_distr_discr_set_pv(distr, pmf, length(pmf))
    par=unur_dgt_new(distr)
    unur_set_urng(par, rng)
    gen=unur_init(par)
    (gen, distr)
end

type DiscretePVUnur
    parameters::Array{Float64,1}
    times::Array{Float64,1}
    te::Float64
    gen::Ptr{UNUR_GEN}
    distr::Ptr{UNUR_DISTR}
end

function DiscretePVUnur(rng::Ptr{UNUR_URNG}, times, pv)
    gen, distr=MakeDiscretePVUnur(rng, pv)
    DiscretePVUnur(Float64[], times, 0.0, gen, distr)
end

Parameters(gu::DiscretePVUnur)=gu.parameters
function Parameters!(gu::DiscretePVUnur, rng)
end

