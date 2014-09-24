include("SmallGraph.jl")

using SmallGraphs

g=DirectedGraph()
add_node(g, 1, {"mark"=>7})
add_edge(g, 3, 5, {"stoch"=>-1})
add_edge(g, 3, 7)

println(g.node[1])
println(g.node[3])
println(g.edge[3][5])
println(neighbors(g, 3))

u=UndirectedGraph()
add_node(u, "s", {"mark"=>99})
add_node(u, "i", {"mark"=>5})
add_node(u, "r", {"mark"=>0})
add_node(u, "infect")
add_node(u, "recover")

add_edge(u, "s", "infect", {"stoch"=>-1})
add_edge(u, "i", "infect", {"stoch"=>-1})
add_edge(u, "infect", "i", {"stoch"=>2})
add_edge(u, "i", "recover", {"stoch"=>-1})
add_edge(u, "recover", "r", {"stoch"=>1})

deps=Dict{Any,Set{Any}}()
deps["s"]=Set{Any}(["infect"])
deps["i"]=Set{Any}(["infect","recover"])
deps["r"]=Set{Any}()
