
using Distributions

import Base: push!, length, pop!, get, print
export GSPNModel, ConstExplicitTransition
export ExplicitGSPN, add_place, add_transition, current_time
export sir_explicit, fire, enabled_transitions, print
export TransitionRoute

include("explicitgraphs.jl")

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


type ExplicitGSPN{PK,TK}
    # Graph with stoichiometry
    gspn::GSGraph{PK,TK}
    # Transition properties starts from index of first transition.
    transition_prop::Array{ConstExplicitTransition,1}
    ExplicitGSPN()=new(GSGraph{PK,TK}(), Array(ConstExplicitTransition,0))
end

function print(io::IO, g::ExplicitGSPN)
    # println(io, "GSPN v ", num_vertices(g.gspn), " e ", num_edges(g.gspn))
    # for (id, name) in enumerate(g.id_to_pt)
    #     assert(g.pt_to_id[name]==id)
    # end
    # print(io, "places ")
    # for ppidx in 1:g.transition_base
    #     print(io, g.id_to_pt[ppidx], ", ")
    # end
    # println(io)

    # println(io, "transitions")
    # assert(length(g.id_to_pt)-g.transition_base==length(g.transition_prop))
    # for ptidx in (g.transition_base+1):length(g.id_to_pt)
    #     name=g.id_to_pt[ptidx]
    #     print(io, name, " ", ptidx, ": ")
    #     transition_places(g, ptidx) do place_id, stoich, name
    #         place=g.id_to_pt[place_id]
    #         print(io, "(", place, ", ", stoich, ", ", name, ") ")
    #     end
    #     println(io)
    # end
end

transitionobj(eg::ExplicitGSPN, id)=eg.transition_prop[id]
transitionids(eg::ExplicitGSPN)=1:length(eg.gspn.transitions)

#### Construction
function add_place(structure::ExplicitGSPN, place)
    add_place!(place, g.gspn)
end

# Will have an error b/c this is used for transitions, too.
function place_to_key(structure::ExplicitGSPN, place)
    place_to_key(place, structure.gspn)
end

function key_to_place(structure::ExplicitGSPN, id::Int)
    key_to_place(id, structure.gspn)
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
    add_transition!(transition_name, structure.gspn)
    push!(structure.transition_prop, transition)

    for (route_idx, route) in enumerate(transition_routes)
        if route.from_place!=nothing
            epf=EdgeProperty(-route.stoichiometry, route_idx, :Any)
            add_edge!(transition_name, route.from_place, epf, structure.gspn)
        end
        if route.to_place!=nothing
            ept=EdgeProperty(route.stoichiometry, route_idx, :Any)
            add_edge!(transition_name, route.to_place, structure.gspn)
        end
    end
    # The dependency graph is a separate entity but
    # saved in the same structure with stoichiometry=0.
    for (place, local_name) in dependencies
        add_dependency!(transition_name, place, structure.gspn)
    end
end

# Iterates over the gspn.
function transition_places(f::Function, eg::ExplicitGSPN, transition_id::Int64)
    for edge in out_edges(transition_id, eg.gspn)
        f(edge.n, edge.property.stoichiometry, edge.property.route)
    end
end

# Return a transition's dependency by its index.
function transition_dependency(eg::ExplicitGSPN, transition_id::Int, idx::Int)
    out_dependency(transition_id, eg.dependency)[idx]
end

function transition_dependency_length(eg::ExplicitGSPN, transition_id::Int)
    out_depdegree(transition_id, eg.dependency)
end

function dependent_transitions(eg::ExplicitGSPN, place_id::Set{Int64})
    # affected=Set{Int64}()
    # for p in place_id
    #     for trans_neighbor in in_neighbors(p, eg.gspn)
    #         push!(affected, trans_neighbor)
    #     end
    # end
    # affected
    foldl(union!, Set{Int}(),
        [in_neighbors(p, eg.gspn) for p in place_id])
end

