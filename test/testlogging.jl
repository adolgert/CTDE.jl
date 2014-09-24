include("semimarkov.jl")

using SemiMarkov

logging_level("debug")

@trace("Trace this", 7)
@debug("Debbb")
@info("infotastic")
@warn("warn me? warn you")
@error("fubar!")
