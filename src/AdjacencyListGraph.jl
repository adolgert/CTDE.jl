module AdjacencyListGraph
using DataStructures

include("container.jl")

import Base: start, done, next, isequal
export AbstractAdjacencyList
# Export the neighbor and vertex because the types are in 
# vertex_descriptor's type sometimes.
export AdjacencyList, AdjacencyListNeighbor, AdjacencyListVertex
export add_vertex!, add_edge!
export graph_property, vertex_property, edge_property
export num_vertices, vertices, num_edges, out_edges
export out_degree, out_neighbors, source, target
export start, done, next, isequal


abstract AbstractAdjacencyList{VP,EP,GP}

immutable type AdjacencyListNeighbor{ALV,EP}
    n::ALV
    edge_property::EP
end

# A Set container should contain one of each vertex, even if properties
# are specified differently.
function isequal(a::AdjacencyListNeighbor, b::AdjacencyListNeighbor)
    isequal(a.n, b.n)
end


# This fools the type system into creating v::Vector{Int} kinds
# of containers. It is dirty pool.
store_construct(x::TypeVar, args...)=x
function store_construct(NC::DataType, VC, EP, self)
    alv=container_key(VC{self})
    NC{AdjacencyListNeighbor{alv,EP}}
end


immutable type AdjacencyListVertex{VP,EP,VC,NC}
    v::store_construct(NC, VC, EP, AdjacencyListVertex{VP,EP,VC,NC})
    vertex_property::VP
    AdjacencyListVertex(vp::VP)=new(
        container_construct(store_construct(
            NC, VC, EP, AdjacencyListVertex{VP,EP,VC,NC})),
        vp)
end


function AdjacencyListVertex{VP}(VC, NC, EP, vp::VP)
    AdjacencyListVertex{VP,EP,VC,NC}(vp)
end


type AdjacencyList{VP,EP,GP,VC,NC} <: AbstractAdjacencyList{VP,EP,GP}
    vertices::cconstruct(VC, AdjacencyListVertex{VP,EP,VC,NC})
    is_directed::Bool
    graph_property::GP
    AdjacencyList(directed::Bool, gp::GP)=new(
        container_construct(cconstruct(
            VC, AdjacencyListVertex{VP,EP,VC,NC})),
        directed,
        gp
        )
end


function AdjacencyList{GP}(VC::Union(DataType,TypeConstructor),  
        NC::Union(DataType,TypeConstructor), VP::Union(DataType,()),
        EP::Union(DataType,()), gp::GP)
    AdjacencyList{VP,EP,GP,VC,NC}(false, gp)
end


function AdjacencyList{GP}(VC::Union(DataType,TypeConstructor),  
        NC::Union(DataType,TypeConstructor), VP::DataType,
        EP::DataType, gp::GP, capacity::Int)
    adj=AdjacencyList{VP,EP,GP,VC,NC}(false, gp)
    for i in 1:capacity
        add_vertex!(adj)
    end
    adj
end


function make_neighbor{VP,EP,GP,VC,NC}(v, ep, g::AdjacencyList{VP,EP,GP,VC,NC})
    container_value(store_construct(
        NC, VC, EP, AdjacencyListVertex{VP,EP,VC,NC}))(
        v, ep)
end


function add_vertex!{VP,EP,GP,VC,NC}(g::AdjacencyList{VP,EP,GP,VC,NC})
    vertex=AdjacencyListVertex(VC, NC, EP, VP())
    container_push(g.vertices, vertex)
end


function add_vertex!{VP,EP,GP,VC,NC}(vp::VP, g::AdjacencyList{VP,EP,GP,VC,NC})
    vertex=AdjacencyListVertex(VC, NC, EP, vp)
    container_push(g.vertices, vertex)
end

# For when vertices are in a dict. k is the key.
function add_vertex!{VP,EP,GP,VC,NC}(k, vp::VP, g::AdjacencyList{VP,EP,GP,VC,NC})
    vertex=AdjacencyListVertex(VC, NC, EP, vp)
    container_push(g.vertices, k, vertex)
end


source(edge, g::AdjacencyList)=edge[1]
function target(edge, g::AdjacencyList)
    edge[2].n
end

# An edge is a tuple of the vertex_descriptor of the source
# and the vertex value of the target.
function edge(u, v, g::AdjacencyList)
    (u, container_get(container_get(g.vertices, u).v, v) )
end

function add_edge!{VP,EP,GP,VC,NC}(u, v, ep::EP, g::AdjacencyList{VP,EP,GP,VC,NC})
    uvertex=container_get(g.vertices, u)
    n=make_neighbor(v, ep, g)
    right=container_push(uvertex.v, n)
    if !g.is_directed
        vvertex=container_get(g.vertices, v)
        left=container_push(vvertex.v, make_neighbor(u, ep, g))
    end
    (u, n) # This is the edge descriptor.
end

# For when edges are in a dict. k is the key to the edge.
function add_edge!{VP,EP,GP,VC,NC}(u, v, k, ep::EP, g::AdjacencyList{VP,EP,GP,VC,NC})
    uvertex=container_get(g.vertices, u)
    n=make_neighbor(v, ep, g)
    right=container_push(uvertex.v, k, n)
    (u, n) # This is the edge descriptor.
end


graph_property(g::AdjacencyList)=g.gp
vertex_property(vertex_descriptor, g::AdjacencyList)=g.vertices[vertex_descriptor].vp
edge_property(edge_descriptor, g::AdjacencyList)=edge_descriptor[2].ep

num_vertices(g::AdjacencyList)=length(g.vertices)

vertices(g::AdjacencyList)=container_iter(g.vertices)

function num_edges{V,E}(g::AdjacencyList{V,E})
    sum([length(x.v) for v in g.vertices])
end

# Needs types
type OutEdgeNeighborIter{N,S}
    neighbor_iter::N
    source::S
end

start(iter::OutEdgeNeighborIter)=start(iter.neighbor_iter)
done(iter::OutEdgeNeighborIter, state)=done(iter.neighbor_iter, state)
function next(iter::OutEdgeNeighborIter, state)
    n=next(iter.neighbor_iter, state)
    (iter.source, n)
end

function out_edges(vertex_descriptor, g::AdjacencyList)
    OutEdgeNeighborIter(container_value_iter(g.vertices[vertex_descriptor].v),
        vertex_descriptor)
end


out_degree(vertex_descriptor, g::AdjacencyList)=
        length(g.vertices[vertex_descriptor].v)


type AdjacencyListNeighborIter{I,V}
    inner::I
    v::V
    AdjacencyListNeighborIter(v)=new(container_iter(v), v)
end

start(iter::AdjacencyListNeighborIter)=start(iter.inner)
function next(iter::AdjacencyListNeighborIter, idx)
    container_get(v, next(iter.inner, idx)).n
end
done(iter::AdjacencyListNeighborIter, idx)=done(iter.inner, idx)

out_neighbors(vertex_descriptor, g::AdjacencyList)=
    AdjacencyListNeighborIter(g.vertices[vertex_descriptor].v)


end # Module
