include("marking.jl")

creator=DefaultTokenCreator{EmptyToken}()
mark=TokenMarking(Dict{Any,EmptyToken}(), creator)
marking=TokenMarking(EmptyToken)
creator=IdTokenCreator()
id_marking=TokenMarking(creator)

println(create(marking.token_creator, 3))
println(create(id_marking.token_creator, 6))
println(create(id_marking.token_creator, 6))

