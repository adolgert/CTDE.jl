module CTDE
    using Logging
    @Logging.configure(level=INFO)
    include("samplesemimarkov.jl")
    include("transitiondistributions.jl")
    include("partialprocess.jl")
    include("category_fsm.jl")
    include("smallgraph.jl")
    have_unuran=false
    try
    	Pkg.installed("UNURAN")
    	have_unuran=true
    end
    if have_unuran
    	include("unudist.jl")
    end
end
