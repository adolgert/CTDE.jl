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
    v::ALV
    edge_property::EP
end

# A Set container should contain one of each vertex, even if properties
# are specified differently.
function isequal(a::AdjacencyListNeighbor, b::AdjacencyListNeighbor)
    isequal(a.v, b.v)
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
        NC::Union(DataType,TypeConstructor), VP::DataType,
        EP::DataType, gp::GP)
    AdjacencyList{VP,EP,GP,VC,NC}(false, gp)
end


function make_neighbor{VP,EP,GP,VC,NC}(v, ep, g::AdjacencyList{VP,EP,GP,VC,NC})
    container_value(store_construct(
        NC, VC, EP, AdjacencyListVertex{VP,EP,VC,NC}))(
        v, ep)
end


function add_vertex!{VP,EP,GP,VC,NC}(vp::VP, g::AdjacencyList{VP,EP,GP,VC,NC})
    vertex=AdjacencyListVertex(VC, NC, EP, vp)
    container_push(g.vertices, vertex)
end


function add_edge!{VP,EP,GP,VC,NC}(u, v, ep::EP, g::AdjacencyList{VP,EP,GP,VC,NC})
    uvertex=container_get(g.vertices, u)
    right=container_push(uvertex.v, make_neighbor(v, ep, g))
    (u, right) # This is the edge descriptor.
end

graph_property(g::AdjacencyList)=g.gp
vertex_property(vertex_descriptor, g::AdjacencyList)=g.vertices[vertex_descriptor].vp
function edge_property(edge_descriptor, g::AdjacencyList)
    container_get(container_get(g.vertices, source(edge_descriptor)).v,
        edge_descriptor[2]).ep
end

num_vertices(g::AdjacencyList)=length(g.vertices)

vertices(g::AdjacencyList)=container_iter(g.vertices)

function num_edges{V,E}(g::AdjacencyList{V,E})
    sum([length(x.v) for v in g.vertices])
end

# Needs types
type OutEdgeNeighborIter
    neighbor_iter
    source
end

start(iter::OutEdgeNeighborIter)=start(iter.neighbor_iter)
done(iter::OutEdgeNeighborIter, state)=done(iter.neighbor_iter, state)
function next(iter::OutEdgeNeighborIter, state)
    n=next(iter.neighbor_iter, state)
    (iter.source, n)
end

function out_edges(vertex_descriptor, g::AdjacencyList)
    OutEdgeNeighborIter(g.vertices[vertex_descriptor].v, vertex_descriptor)
end


out_degree(vertex_descriptor, g::AdjacencyList)=
        length(g.vertices[vertex_descriptor].v)


type AdjacencyListNeighborIter
    inner
    v
    AdjacencyListNeighborIter(v)=new(iter(v), v)
end

start(iter::AdjacencyListNeighborIter)=start(iter.inner)
function next(iter::AdjacencyListNeighborIter, idx)
    container_get(v, next(iter.inner, idx)).v
end
done(iter::AdjacencyListNeighborIter, idx)=done(iter.inner, idx)

out_neighbors(vertex_descriptor, g::AdjacencyList)=
    AdjacencyListNeighborIter(g.vertices[vertex_descriptor].v)

source(edge, g::AdjacencyList)=edge[1]
function target(edge, g::AdjacencyList)
    container_get(container_get(g.vertices, edge[1]), edge[2]).v
end

end # Module
