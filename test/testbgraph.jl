using DataStructures
push!(LOAD_PATH, "../../src")
using CTDE.AdjacencyListGraph
nei=AdjacencyListGraph.AdjacencyListNeighbor(3, 17)

adj=AdjacencyListGraph.AdjacencyList(Vector, Deque, Int, Float64, "my graph")
v1=AdjacencyListGraph.add_vertex!(24, adj)
v2=AdjacencyListGraph.add_vertex!(25, adj)
e1=AdjacencyListGraph.add_edge!(v1, v2, 3.7, adj)
println(e1)
