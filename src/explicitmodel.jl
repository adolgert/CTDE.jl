
using Distributions
using Graphs

import Base: push!, length, pop!, get, print
export GSPNModel, ConstExplicitTransition
export ExplicitGSPN, add_place, add_transition, current_time
export sir_explicit, fire, enabled_transitions, print
export TransitionRoute


abstract ExplicitTransition
# This is const because it has no internal state
# besides the two functions.
type ConstExplicitTransition <: ExplicitTransition
    distribution::Function
    fire::Function
    ConstExplicitTransition(d,f)=new(d, f)
end

ConstExplicitTransition(dist)=ConstExplicitTransition(dist, (x...)->x)


# Calculates the distribution, given marking at the places
# upon which this transition depends.
# Call through method in case a transition saves state.
# A coalesced transition would need internal state.
function distribution(et::ExplicitTransition, local_marking, user, when)
    et.distribution(local_marking, user, when)
end

function fire(et::ExplicitTransition, tokendict, user)
    et.fire(tokendict, user)
end


immutable type EdgeProperty
    stoichiometry::Int
    route::Int
    selection::Symbol
end

# Says some tokens go from this place to that one.
# They can be modified in-between.
# Used during specification of the transition.
type TransitionRoute
  from_place
  to_place
  stoichiometry::Int
  selection::Symbol
  # Whether this route can create or destroy tokens.
  create_or_destroy::Int
  TransitionRoute(f, t, s)=new(f, t, s, :Any, 0)
end

typealias DependencyGraph GenericAdjacencyList{Int64,Array{Int64,1},Array{Array{Int64,1},1}}

type ExplicitGSPN
    # Graph with stoichiometry
    gspn::IncidenceList{Int64,Edge{Int64}}
    dependency::DependencyGraph
    pt_to_id::Dict{Any,Int64}
    id_to_pt::Array{Any,1}
    edge_prop::Array{EdgeProperty,1}
    # Transition properties starts from index of first transition.
    transition_prop::Array{ConstExplicitTransition,1}
    # Offset of transitions in the gspn.
    transition_base::Int64
end

function ExplicitGSPN()
    gspn=inclist(Int64, is_directed=false)
    dependency=adjlist(Int, is_directed=true)
    ExplicitGSPN(gspn, dependency, Dict{Any,Int64}(),Array(Any,0),
        Array(EdgeProperty,0), Array(ConstExplicitTransition,0), -1)
end

function print(io::IO, g::ExplicitGSPN)
    println(io, "GSPN v ", num_vertices(g.gspn), " e ", num_edges(g.gspn))
    for (id, name) in enumerate(g.id_to_pt)
        assert(g.pt_to_id[name]==id)
    end
    print(io, "places ")
    for ppidx in 1:g.transition_base
        print(io, g.id_to_pt[ppidx], ", ")
    end
    println(io)

    println(io, "transitions")
    assert(length(g.id_to_pt)-g.transition_base==length(g.transition_prop))
    for ptidx in (g.transition_base+1):length(g.id_to_pt)
        name=g.id_to_pt[ptidx]
        print(io, name, " ", ptidx, ": ")
        transition_places(g, ptidx) do place_id, stoich, name
            place=g.id_to_pt[place_id]
            print(io, "(", place, ", ", stoich, ", ", name, ") ")
        end
        println(io)
    end
end

transitionobj(eg::ExplicitGSPN, id)=eg.transition_prop[id-eg.transition_base]
transitionids(eg::ExplicitGSPN)=(eg.transition_base+1):length(eg.id_to_pt)

#### Construction
function add_place(structure::ExplicitGSPN, place)
    if structure.transition_base>=0
        error("Please add places before adding any transitions.")
    end
    push!(structure.id_to_pt, place)
    id=length(structure.id_to_pt)
    structure.pt_to_id[place]=id
    add_vertex!(structure.gspn, id)
    add_vertex!(structure.dependency, id)
end

function place_to_key(structure::ExplicitGSPN, place)
    structure.pt_to_id[place]
end

function key_to_place(structure::ExplicitGSPN, id::Int)
    structure.id_to_pt[id]
end

# id=transition key
# transition=ExplicitTransition object
# Route - See TransitionRoute type. Tokens move along routes betwen places.
# dependencies=iterable of (place, local name)
# The local name for stoichiometry is used in the transition.fire
# method to find input tokens, no matter which place is hooked
# to the transition.
# The local name for dependencies is used in the transition.distribution
# method.
function add_transition(structure::ExplicitGSPN, transition_name, transition,
        transition_routes, dependencies)
    if structure.transition_base<0
        structure.transition_base=length(structure.id_to_pt)
    end
    push!(structure.id_to_pt, transition_name)
    id=length(structure.id_to_pt)
    structure.pt_to_id[transition_name]=id
    add_vertex!(structure.gspn, id)
    add_vertex!(structure.dependency, id)
    push!(structure.transition_prop, transition)

    for (route_idx, route) in enumerate(transition_routes)
        if route.from_place!=nothing
            from_place_id=structure.pt_to_id[route.from_place]
            add_edge!(structure.gspn, id, from_place_id)
            epf=EdgeProperty(-route.stoichiometry, route_idx, :Any)
            push!(structure.edge_prop, epf)
        end
        if route.to_place!=nothing
            to_place_id=structure.pt_to_id[route.to_place]
            add_edge!(structure.gspn, id, to_place_id)
            ept=EdgeProperty(route.stoichiometry, route_idx, :Any)
            push!(structure.edge_prop, ept)
        end
    end
    # The dependency graph is a separate entity but
    # saved in the same structure with stoichiometry=0.
    for (place, local_name) in dependencies
        place_id=structure.pt_to_id[place]
        add_edge!(structure.dependency, id, place_id)
    end
end

# Iterates over the gspn.
function transition_places(f::Function, eg::ExplicitGSPN, transition_id::Int64)
    for edge in out_edges(transition_id, eg.gspn)
        ep=eg.edge_prop[edge_index(edge)]
        f(target(edge, eg.gspn), ep.stoichiometry, ep.route)
    end
end

# Return a transition's dependency by its index.
function transition_dependency(eg::ExplicitGSPN, transition_id::Int, idx::Int)
    out_neighbors(transition_id, eg.dependency)[idx]
end

function transition_dependency_length(eg::ExplicitGSPN, transition_id::Int)
    out_degree(transition_id, eg.dependency)
end

function dependent_transitions(eg::ExplicitGSPN, place_id::Set{Int64})
    affected=Set{Int64}()
    for p in place_id
        for oe in out_edges(p, eg.gspn)
            push!(affected, target(oe, eg.gspn))
        end
    end
    affected
end

