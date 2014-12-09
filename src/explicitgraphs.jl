
immutable type EdgeProperty
    stoichiometry::Int
    route::Int
    selection::Symbol
end

immutable type GSNeighborPlace
    n::Int
    property::EdgeProperty
end


type PlaceVertex{PlaceKey}
    r::Vector{Int} # In edges from transitions.
    key::PlaceKey
    PlaceVertex(key::PlaceKey)=new(Array(Int,0), key)
end

type TransitionVertex{TransitionKey}
    v::Vector{GSNeighborPlace} # Out edges to places.
    dep::Vector{Int} # Out edges for dependency.
    key::TransitionKey
    TransitionVertex(key::TransitionKey)=new(
        Array(GSNeighborPlace, 0),
        Array(Int, 0),
        key)
end


type GSGraph{PlaceKey,TransitionKey}
    places::Vector{PlaceVertex{PlaceKey}}
    transitions::Vector{TransitionVertex{TransitionKey}}
    place_key_to_int::Dict{PlaceKey,Int}
    transition_key_to_int::Dict{TransitionKey,Int}
    GSGraph()=new(Array(PlaceVertex{PlaceKey}, 0),
        Array(TransitionVertex{TransitionKey}, 0),
        Dict{PlaceKey,Int}(),
        Dict{TransitionKey,Int}())
end


function add_place!{PK,TK}(place_key::PK,  g::GSGraph{PK,TK})
    push!(g.places, PlaceVertex{PK}(place_key))
    place_key_to_int[place_key]=length(g.places)
end

function add_transition!{PK,TK}(transition_key::TK,  g::GSGraph{PK,TK})
    push!(g.transitions, TransitionVertex{TK}(transition_key))
    transition_key_to_int[transition_key]=length(g.transitions)
    length(g.transitions)
end

function add_edge!{PK,TK}(transition::TK, place::PK,
        ep::EdgeProperty, g::GSGraph{PK,TK})
    tid=transition_key_to_int[transition]
    pid=place_key_to_int[place]

    push!(g.transitions[tid].v, GSNeighborPlace(pid, ep))
    push!(g.places[pid].r, tid)
end

place_to_key{PK,TK}(place::PK, g::GSGraph{PK,TK})=g.place_key_to_int[place]::Int
key_to_place{PK,TK}(id::Int, g::GSGraph{PK,TK})=g.places[id].key::PK
transition_to_key{PK,TK}(transition::TK, g::GSGraph{PK,TK})=
        g.transition_key_to_int[transition]::Int
transition_to_place{PK,TK}(id::Int, g::GSGraph{PK,TK})=g.transitions[id].key::TK

function add_dependency!{PK,TK}(transition::TK, place::PK, g::GSGraph{PK,TK})
    tid=transition_key_to_int[transition]
    pid=place_key_to_int[place]

    push!(g.transitions[tid].dep, pid)
end

function out_edges{PK,TK}(tid::Int, g::GSGraph{PK,TK})
    g.transitions[tid].v
end

function in_neighbors{PK,TK}(pid::Int, g::GSGraph{PK,TK})
    g.places[pid].r
end

function out_dependency{PK,TK}(tid::Int, g::GSGraph{PK,TK})
    g.transitions[tid].dep
end

out_depdegree{PK,TK}(tid::Int, g::GSGraph{PK,TK})=length(g.transitions[tid].dep)

