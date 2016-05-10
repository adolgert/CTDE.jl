# Average and variance of ribosome leaving mRNA as a
# function of time, in absence of mRNA decay.

# Distribution of number of proteins produced

# Probability density of the waiting time between the production
# of E and E-1, with inset of E_0.

using Gadfly

include("protein_copy_number_distribution.jl")


function DecayRuns()
    #rng=MersenneTwister(333333)
    seed=UInt64[111, 222, 333, 444, 555, 666]
    RngStream_SetPackageSeed(seed)
    rng=unur_urng_rngstream_new("urng-1")
    run_cnt=1000
    parameters=Dict(
        :n => 0.5,
        :k => 0.5,
        :alpha => 1/45,
        :beta => 2/15,
        :lambda => 1/1350,
        :L => 815,
        :l => 12
        )
    process, state=MakeProcess(parameters, rng)
    start_time=now()
    obs=Observations()
    for run_idx = 1:run_cnt
        Reset!(state)
        sampler=NextReactionHazards()

        push!(obs.row, length(obs.leave)+1)
        RunSimulation(process, sampler, Observer(obs), rng)
        if now()-start_time > Base.Dates.Minute(30)
            Write(obs, "$(run_cnt)-$(run_idx)")
            Clear!(obs)
            start_time=now()
        end
    end
    Write(obs, "$(run_cnt)")
end


function LongRuns()
    #rng=MersenneTwister(333333)
    seed=UInt64[111, 222, 333, 444, 555, 666]
    RngStream_SetPackageSeed(seed)
    rng=unur_urng_rngstream_new("urng-1")
    run_cnt=10
    parameters=Dict(
        :n => 0.5,
        :k => 0.5,
        :alpha => 1/45,
        :beta => 2/15,
        :lambda => 1/1350,
        :L => 815,
        :l => 12
        )
    process, state=MakeProcess(parameters, rng, false)
    start_time=now()
    obs=Observations()
    for run_idx = 1:run_cnt
        Reset!(state)
        sampler=NextReactionHazards()

        push!(obs.row, length(obs.leave)+1)
        RunSimulation(process, sampler, TimedObserver(obs), rng)
        if now()-start_time > Base.Dates.Minute(30)
            Write(obs, "long$(run_cnt)-$(run_idx)")
            Clear!(obs)
            start_time=now()
        end
    end
    Write(obs, "long$(run_cnt)")
end


function Plot()
	obs=FromFile("1000")
    (hist, proteins_produced)=DistributionOfProteinsProduced(obs)
    # We should lose the value at 0 to match the graph.
    if length(proteins_produced)>51
    	proteins_produced=proteins_produced[2:51]
    else
    	proteins_produced=proteins_produced[2:length(proteins_produced)]
    end
    distplot=plot(x=Array(1:(length(proteins_produced))), y=proteins_produced,
    		Geom.point, Scale.y_log,
    		Guide.Title("Distribution of Number Produced"),
	    	Guide.XLabel("Count"), Guide.YLabel("Log of Fraction"))
    draw(PDF("distribution.pdf", 4inch, 3inch), distplot)
end

DecayRuns()
#LongRuns()
