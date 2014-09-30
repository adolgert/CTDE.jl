module SemiMarkov
    using Logging
    @Logging.configure(level=INFO)
    include("smallgraph.jl")
    include("samplesemimarkov.jl")
    include("marking.jl")
    include("explicitmodel.jl")
    include("transitiondistributions.jl")
end
