using UNURAN

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
	println("Truncate to $left")
	unur_distr_cont_set_domain(distr, left, Inf)
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
	println("SampleThenQuantile sample $v survival $S")
	(v, S)
end

# This one samples by inversion.
function SampleByInversion(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR},
		rng::Ptr{UNUR_URNG})
	U=unur_urng_sample(rng)
	v=unur_quantile(gen, U)
	(v, 1-U)
end

function ConditionalSurvival(distr::Ptr{UNUR_DISTR}, a, b)
	cdf=unur_distr_cont_get_cdf(distr)
	println("Have the cdf")
	if a>eps(Float64)
    	last=ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}),
    			a, distr)
    else
    	last=0.0
    end
	current=ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}),
			b, distr)
	println("ConditionalSurvival a $a b $b last $last current $current")
	(1.0-last)/(1.0-current)
end

# We are given the conditional survival for a distribution
# from some time (called shift) to a future time.
# 1. Find the survival of the distribution to that shifted time.
# 2. Use that and the conditonal survival to invert the distribution.
function ShiftedQuantile(gen::Ptr{UNUR_GEN}, distr::Ptr{UNUR_DISTR},
		survival, shift)
	cdf=unur_distr_cont_get_cdf(distr)
  	s_prime=ccall(cdf, Cdouble, (Cdouble, Ptr{UNUR_DISTR}), shift, distr)
  	unur_quantile(gen, 1.0-survival*s_prime)
end

function randexp(rng::Ptr{UNUR_URNG})
	unur_sample_cont(g_Unit_Exponential::Ptr{UNUR_GEN})
end



export GammaUnur

# Gamma Distribution, implemented with UNURAND
# http://statmath.wu.ac.at/unuran/doc/unuran.html#gamma
type GammaUnur <: TransitionDistribution
	parameters::Array{Float64,1}
	te::Float64
	gen::Ptr{UNUR_GEN}
	distr::Ptr{UNUR_DISTR}
end


function GammaUnur(rng::Ptr{UNUR_URNG}, alpha, beta)
	println("GammaUnur begin")
	distr=unur_distr_gamma([alpha, beta], 2)
	par=unur_ninv_new(distr)
	unur_set_urng(par, rng)
	gen=unur_init(par)
	gu=GammaUnur([alpha, beta], 0.0, gen, distr)
	println("GammaUnur end")
	gu
end


Parameters(gu::GammaUnur)=gu.parameters
function Parameters!(gu::GammaUnur, rng, alpha, beta)
	println("GammaUnur::Parameters! begin")
	gu.parameters[1]=alpha
	gu.parameters[2]=beta
	if gu.gen!=nothing
		unur_free(gu.gen)
		unur_distr_free(gu.distr)
	end
	println("GammaUnur::Create a new one $alpha $beta")
	gu.distr=unur_distr_gamma(Float64[alpha, beta], 2)
	out_params=Ref{Ptr{Cdouble}}(0)
	nparams=unur_distr_cont_get_pdfparams(gu.distr, out_params)
	println("outparams $(out_params[])")
	arr=pointer_to_array(out_params[], nparams)
	println("GammaUnur params: $arr")
	println("GammaUnur::Create parameters $nparams")
	par=unur_auto_new(gu.distr)
	println("GammaUnur::Set urng")
	unur_set_urng(par, rng)
	println("GammaUnur::Init")
	gu.gen=unur_init(par)
	println("GammaUnur::Parameters! end")
end


function EnablingTime!(gu::GammaUnur, t::Float64)
	println("enabling time $t")
    gu.te=t
end


function Sample(dist::GammaUnur, now::Float64, rng)
    d=now-dist.te
    if abs(d)<eps(1.0)
    	println("Sample: te $(gu.te) now $now")
    	v=dist.te+unur_sample_cont(dist.gen)
    else
    	println("Sample shifted te $(dist.te) now $now")
    	TruncateDomain(dist.gen, d)
    	v=dist.te+unur_sample_cont(dist.gen)
    end
    v
end


function MeasuredSample(gu::GammaUnur, now::Float64, rng)
    d=now-gu.te
    if abs(d)<eps(1.0)
    	(v, S)=SampleThenQuantile(gu.gen, gu.distr)
    else
    	println("MeasuredSample: te $(gu.te) now $now")
    	TruncateDomain(gu.gen, d)
    	(v, S)=SampleThenQuantile(gu.gen, gu.distr)
    end
	(gu.te+v, S)
end


# xa is the conditonal survival
function ConsumeSample(gu::GammaUnur, xa, start, finish)
    xa=(xa<0) ? 1 : xa
    xa*ConditionalSurvival(gu.distr, start-gu.te, finish-gu.te)
end


function Putative(gu::GammaUnur, when, interval, consumed_interval)
	# If you want unur_quantile to work, you have to use
	# HINV, NINV, PINV, CSTD, or DGT. CSTD requires that
	# the generator have an inversion method.
	S=interval/consumed_interval
	shift=when-gu.te
	println("Putative S $S shift $shift")
	gu.te+ShiftedQuantile(gu.gen, gu.distr, S, shift)
end
