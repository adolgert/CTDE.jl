
macro warn(msg...)
    apply(println, eval(msg))
end

macro info(msg...)
    apply(println, eval(msg))
end

macro debug(msg...)
    #apply(println, eval(msg))
end

macro trace(msg...)
    #apply(println, eval(msg))
end

# function logging_level(level)
#     levels=["error", "warn", "info", "debug", "trace"]
#     which=findfirst(levels, level)
#     if which==0
#         error("tracing.logging_level ", level, " not among ", levels)
#     end

#     println("logging_level")
#     for lidx=1:length(levels)
#         name=levels[lidx]
#         println("defining ", name)
#         if lidx<=which
#             @eval begin
#                 macro ($name)(args...)
#                     apply(println, eval(args))
#                 end
#             end
#         else
#             @eval begin
#                 macro ($name)(args...)
#                 end
#             end
#         end
#     end
#     @info("Logging level is ", levels[which])
# end
