
using Distributions
using Graphs

import Base: push!, length, pop!, get, print
export ExplicitGSPNModel, ConstExplicitTransition
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
function distribution(et::ExplicitTransition, local_marking, when)
    et.distribution(local_marking, when)
end

function fire(et::ExplicitTransition, tokendict)
    et.fire(tokendict)
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

# Iterates over the dependency graph.
function transition_dependencies(f::Function, eg::ExplicitGSPN, transition_id::Int64,
        lm::Dict{Int,Any})
    for (idx, place) in enumerate(out_neighbors(transition_id, eg.dependency))
        f(place, idx, lm)
    end
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

abstract Model

# A model has both state and structure.
# This is what is presented to the propagator.
type ExplicitGSPNModel{T} <: Model
    structure::ExplicitGSPN
    state::T
end

function ExplicitGSPNModel(state)
    ExplicitGSPNModel(ExplicitGSPN(), state)
end


function ExplicitGSPNModel()
    ExplicitGSPNModel(ExplicitGSPN(), nothing)
end

function current_time(model::ExplicitGSPNModel)
    model.state.current_time
end

function current_time!(model::ExplicitGSPNModel, when::Float64)
    @debug("ExplicitModel.current_time! was ",
        model.state.current_time, " is ", when)
    model.state.current_time=when
end

# This exists to translate between user-defined place names
# and the internal place names. It is a facade.
function add_tokens(model::ExplicitGSPNModel, place::Any, n::Int64)
    place_id=model.structure.pt_to_id[place]
    add_tokens(model.state.marking, place_id, n)
end

function stoichiometry_satisfied(model::ExplicitGSPNModel, transition_id::Int64)
    #@debug("stoichiometry_satisfied enter ", transition_id)
    # Policy: Whether the output places must have no tokens.
    satisfied::Bool=true
    transition_places(model.structure,
            transition_id) do place::Int64, stoich::Int64, route::Int
        if stoich<0
            token_cnt=length(model.state.marking, place)::Int64
            @debug("stoich ", transition_id, " (", place, ", ", stoich, ", ", token_cnt, ")")
            if (stoich+token_cnt)<0
                satisfied=false
            end
        end
    end
    satisfied
end


# This is a way to give a transition access to its dependencies.
# We don't copy them somewhere, just access them when they are called.
type LocalMarking
    dependency::DependencyGraph
    transition_id::Int
    marking
end

function get(lm::LocalMarking, index::Int)
    lm.marking[out_neighbors(lm.transition_id, lm.dependency)[index]]
end

length(lm::LocalMarking)=length(out_neighbors(lm.transition_id, lm.dependency))


# The correct enabling time to send is the current time if we are asking whether
# this transition is newly-enabled. If we are just getting the distribution of
# an enabled transition, then send its already-known enabling time.
function transition_distribution(model::ExplicitGSPNModel, transition_id, enabling)
    local_marking=LocalMarking(model.structure.dependency, transition_id,
        model.state.marking)
    transition=transitionobj(model.structure, transition_id)
    dist, invariant=distribution(transition, local_marking, enabling)
end

function examine_transition(model::ExplicitGSPNModel, transition_id,
        enable::Function, disable::Function)
    #@debug("examine_transition enter ", transition_id)
    stoichiometrically_allowed=stoichiometry_satisfied(model, transition_id)
    @debug("examine_transition allowed ", transition_id, " ",
            stoichiometrically_allowed)
    if !stoichiometrically_allowed
        if haskey(model.state.enabling, transition_id)
            pop!(model.state.enabling, transition_id)
            disable(transition_id, current_time(model))
        else
            # Wasn't enabled. Still not enabled.
        end
        return
    end

    dist, invariant=transition_distribution(model, transition_id, current_time(model))
    nonzero_hazard=(dist!=nothing)
    was_enabled=haskey(model.state.enabling, transition_id)
    if was_enabled
        enabling_record=model.state.enabling[transition_id]
        if nonzero_hazard
            # If there's no invariant, then it's enabled whenever it's on.
            if enabling_record.invariant!=invariant
                enabling_record.time=current_time(model)
                enabling_record.invariant=invariant
                enable(transition_id, dist, current_time(model))
            else
                # Was enabled. Still enabled, same value.
            end
        else
            pop!(model.state.enabling, transition_id)
            disable(transition_id, current_time(model))
        end
    else
        if nonzero_hazard
            model.state.enabling[transition_id]=EnablingRecord(
                    current_time(model),invariant)
            enable(transition_id, dist, current_time(model))
        else
            # Was disabled. Still disabled.
        end
    end
    #@debug("examine_transition exit")
end

init(model::ExplicitGSPNModel)=enable_transitions(model)

# Enables transitions consistent with the current marking.
function enable_transitions(model::ExplicitGSPNModel)
    for transition_id in transitionids(model.structure)
        examine_transition(model, transition_id, (x...)->nothing, (x...)->nothing)
    end
    @debug("enable_transitions enabled ", model.state.enabling)
end


function enabled_transitions(report::Function, model::ExplicitGSPNModel)
    for (id, enabling) in model.state.enabling
        # This isn't asking who becomes enabled but who is already enabled.
        # So use the already-given enabling time.
        distribution, invariant=transition_distribution(model, id, enabling.time)
        @assert(distribution!=nothing)
        report(id, distribution, current_time(model))
    end
    @debug("enabled_transitions enabled ", model.state.enabling)
end


function transition_action(model::ExplicitGSPNModel, transition_id)
    affected_places=Set{Int64}()
    tokens=Dict{Int,Any}() # Map from the edge to the tokens at a place on the edge.
    seen_routes=Set{Int}()
    transition_places(model.structure, transition_id) do place, stoich, route
        if stoich<0
            tokens[route]=empty_container(model.state.marking)
            take!(model.state.marking, place, tokens[route],
                    -stoich)
            push!(affected_places, place)
            push!(seen_routes, route)
        # If a route has no from_place, we still need to give it a spot.
        elseif (stoich>0) && !(route in seen_routes)
            tokens[route]=empty_container(model.state.marking)
        end
    end
    
    @debug("SimpleEGSPN.fire in-tokens ", tokens)
    fire(transitionobj(model.structure, transition_id), tokens)
    # fire_function=transitionobj(model.structure, transition_id).fire
    # fire_function(tokens)

    @debug("SimpleEGSPN.fire out-tokens ", tokens)
    transition_places(model.structure, transition_id) do place, stoich, route
        if stoich>0
            # Tokens can be created here.
            fill!(model.state.marking, place, tokens[route], stoich)
            push!(affected_places, place)
        end
    end
    # Policy: What to do with too many tokens.
    affected_places
end

function fire(model::ExplicitGSPNModel, id_time::NRTransition,
        enable::Function, disable::Function, rng)
    @debug("ExplicitGSPNModel.fire enter ", model.state.last_fired, " idt ",
            id_time)
    transition_id=id_time.key
    affected_places=transition_action(model, transition_id)

    mark_disabled!(model.state, transition_id)
    model.state.last_fired=model.structure.id_to_pt[transition_id]
    current_time!(model, id_time.time)

    examine_transitions=dependent_transitions(model.structure, affected_places)
    # @debug("ExplicitGSPNModel.fire affected ", affected_places,
    #         " exam_trans ", examine_transitions)
    for t in examine_transitions
        examine_transition(model, t, enable, disable)
    end
end

function fire(model::ExplicitGSPNModel, id_time::NRTransition, rng)
    fire(model, id_time, (x...)->nothing, (x...)->nothing, rng)
end

