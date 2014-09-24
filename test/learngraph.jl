using Graphs

v=[1, 2, 3, 5, 7]
verts=ExVertex[ExVertex(x,"") for x in v]
g=inclist(ExVertex, ExEdge{ExVertex})
labels=["s","i","r"]
vertices=Array(ExVertex,3)
for i in 1:3
	vertices[i]=ExVertex(i,labels[i])
	add_vertex!(g, vertices[i])
end
add_vertex!(g, ExVertex(7, "infect"))
add_vertex!(g, ExVertex(9, "recover"))

edge_idx=1
si=ExEdge(edge_idx, vertices[1], ExVertex(2, "infect"), (UTF8String=>Any)["s"=>-1])
add_edge!(g, si)
ii=ExEdge(edge_idx, ExVertex(7, "infect"), vertices[3], (UTF8String=>Any)["s"=>-1])
add_edge!(g, si)

for v in g.vertices
	println(v)
end

vi=vertex_index(ExVertex(7,"blah"), g)
println(typeof(vi))
println(vi)
println(out_edges(ExVertex(1,"duh"),g))
println(out_edges(ExVertex(7,"duh"),g))

# The vertex index is used both as a unique label and as the offset
# into the string of indices, so I have to set it at the moment it goes in.
# Then, with luck, it will coincide.
# In conclusion, Graphs implementation is so confused that it's worthless.
