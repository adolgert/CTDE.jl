module CTDE
    using Logging
    @Logging.configure(level=INFO)
    include("prefixsearch.jl")
    include("sample/direct.jl")
    # include("samplesemimarkov.jl")
    # include("transitiondistributions.jl")
    # include("partialprocess.jl")
    # include("category_fsm.jl")
    # include("smallgraph.jl")
end
