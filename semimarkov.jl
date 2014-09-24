module SemiMarkov
    include("tracing.jl")
    include("smallgraph.jl")
    #include("SimpleExponentialGSPN.jl")
    include("samplesemimarkov.jl")
    include("marking.jl")
    include("explicitmodel.jl")
    include("transitiondistributions.jl")
end
