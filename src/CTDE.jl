module CTDE
    using Logging
    @Logging.configure(level=INFO)
    include("samplesemimarkov.jl")
    include("transitiondistributions.jl")
    include("partialprocess.jl")
    include("category_fsm.jl")
end
