# This is a set of functions to test the TransitionDistribution
# types. Making your own distributions is a bad idea, so let's
# ensure it isn't a damning idea.
#
# Each distribution has three kinds of parameters
# 1. Basic parameters, such as lambda, k for Weibull.
#    There may be a map from the TransitionDistribution parameters
#    to the parameters used by Julia distributions, for instance
#    a rate for a scale=1/rate.
# 2. An optional shift of a distribution to the right,
#    so the hazard is zero for some time before the start.
# 3. An enabling time, which places the distribution on an absolute
#    time scale.
#
# Our tests will include
# 1. Plot the CDF and a kernel density estimation of the PDF
#    for a few paradigmatic values, for instance those in
#    the Wikipedia or Mathematica web pages for comparison.
# 2. A set of statistics, such as mean, variance, skew, and kurtosis
#    are often known values from the parameters. Sometimes it's hard
#    to find variance when these are empirical from a distribution,
#    but it's a start.
# 3. A known good distribution with which to compare using
#    a univariate continuous distribution test, such as
#    Kolmogorov-Smirnov, Cramer-von Mises,
#    or Anderson-Darling. This has a few parts:
#
#    a) Test distribution starting at Te=0, shift=0.
#    b) Sample that distribution in a series of steps (as
#       used in the Next Reaction Method).
#    c) Translate the sample with a Te and compare.
#    d) Translate the sample with a right shift, if appropriate.

using Logging
using Distributions
using Gadfly
using KernelDensity
@Logging.configure(level=INFO)

include("transitiondistributions.jl")

# Define a set of sampling methods for a distribution
function SingleSample(d::TransitionDistribution, N::Int, rng)
	ed=EmpiricalDistribution()
	for i = 1:N
		push!(ed, MeasuredSample(d, 0.0, rng)[1])
	end
	ed
end


# This emulates the set of steps to sample a Next Reaction
# algorithm.
function MultiSample(d::TransitionDistribution, N::Int, rng)
	ed=EmpiricalDistribution()
	for i = 1:N
		consumed=-1
		sample, quantile=MeasuredSample(d, 0.0, rng)
		sample_cnt=5
		for consume_idx=1:sample_cnt
			t=(consume_idx*0.9*sample)/sample_cnt
			tm1=((consume_idx-1)*0.9*sample)/sample_cnt
			quantile=ConsumeSample(d, consumed, tm1, t)
		end
		push!(ed, Putative(d, 0.9*sample, quantile, consumed))
	end
	ed
end


function GraphAt(d::TransitionDistribution, params, name, rng)
	sample_cnt=10000
	xlimit=2.5

	colors=[colorant"blue", colorant"orange", colorant"magenta",
		colorant"green"]
	kde_dist=Distributions.Uniform(-0.05, 0.05)
	pdflayers=[]
	cdflayers=[]


	for idx = 1:length(params)
		Parameters!(d, params[idx]...)
		ed=SingleSample(d, sample_cnt, rng)
		build!(ed)
		xcnt=searchsortedlast(ed.samples, xlimit)
		kernel=kde(ed.samples, boundary=(0.0, 10.0), kde_dist)
		pdfx=Array{Int,1}(0:500)
		pdfx=pdfx/500.0
		theme=Theme(default_color=colors[idx])
		push!(pdflayers, layer(x=pdfx, y=KernelDensity.pdf(kernel, pdfx),
				Geom.line, theme))
		push!(cdflayers, layer(x=CDFX(ed)[1:xcnt], y=CDFY(ed)[1:xcnt],
				Geom.line, theme))
	end
	pdf=plot(pdflayers..., Guide.XLabel("Time"),
			Scale.x_continuous(minvalue=0.0, maxvalue=xlimit),
			Guide.YLabel("Density"),
			Guide.Title("Density Function $name"))
	cdf=plot(cdflayers..., Guide.XLabel("Time"),
			Scale.x_continuous(minvalue=0.0, maxvalue=xlimit),
			Guide.YLabel("PDF"),
			Guide.Title("Cumulative Distribution $name"))
	draw(PDF("$name.pdf", 6inch, 6inch), vstack(pdf, cdf))
end


"""
This prints a set of statistics to check whether a distribution
converges to the expected value.
Each TransitionDistribution defines a set of statistics
that can be used as tests. This generates, say, 10 samples
from the distribution at a particular size, say 10000, and
then asks what the mean and standard deviation of those
10 statistics are. It prints that list so you can see
it converge, or not.
"""
function ProgressiveStatistics(d::TransitionDistribution,
		sampling_function::Function, rng)
	trial_cnt=10
	sample_size=[100, 1000, 10000]
	allstats=Dict{AbstractString,Array}()
	allexpect=Dict{Any,Any}()
	for sample_idx = 1:length(sample_size)
		stats=Dict{AbstractString,Array}()
		for trial_idx = 1:trial_cnt
			ed=sampling_function(d, sample_size[sample_idx], rng)
			onestat=TestStatistics(d, ed)
			for k in keys(onestat)
				allexpect[k]=onestat[k][1]
				if haskey(stats, k)
					push!(stats[k], onestat[k][2])
				else
					stats[k]=[onestat[k][2]]
				end
			end
		end

		# Now add the means and variances to the list
		# for all sample sizes.
		for k in keys(stats)
			m=Base.mean(stats[k])
			std_dev=Base.std(stats[k])
			if haskey(allstats, k)
				push!(allstats[k], [m, std_dev])
			else
				allstats[k]=Any[[m, std_dev]]
			end
		end
	end

	stat_names=AbstractString[]
	for k in keys(allstats)
		push!(stat_names, k)
	end
	sort!(stat_names)
	for name in stat_names
		print("$name=$(allexpect[name])\n")
		for samp_idx = 1:length(sample_size)
			(dmean, dvar)=allstats[name][samp_idx]
			diff=abs(allexpect[name]-dmean)/allexpect[name]
			print("  $(sample_size[samp_idx])\t$dmean\t$diff\t$dvar\n")
		end
	end
end


function Test()
	rng=MersenneTwister(333334)
	weibull=TransitionWeibull(1.0, 2.0)
	params=Any[[1, 0.5], [1, 1], [1, 1.5], [1, 5]]
	GraphAt(weibull, params, "Weibull", rng)

	ProgressiveStatistics(weibull, SingleSample, rng)
end


Test()
