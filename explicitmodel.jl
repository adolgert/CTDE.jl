
using Distributions
using Graphs

import Base: push!, length, pop!, get
export ExplicitGSPNModel, ConstExplicitTransition
export ExplicitGSPN, add_place, add_transition, current_time
export sir_explicit, fire, enabled_transitions


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
    stoichiometry::Int64
    local_name::ASCIIString
end

type ExplicitGSPN
    # Graph with stoichiometry
    gspn::IncidenceList{Int64,Edge{Int64}}
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
    ExplicitGSPN(gspn, Dict{Any,Int64}(),Array(Any,0),
        Array(EdgeProperty,0), Array(ConstExplicitTransition,0), -1)
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
end

# id=transition key
# transition=ExplicitTransition object
# stoichiometry=iterable of (place, +-1, local name)
# dependencies=iterable of (place, local name)
# The local name for stoichiometry is used in the transition.fire
# method to find input tokens, no matter which place is hooked
# to the transition.
# The local name for dependencies is used in the transition.distribution
# method.
function add_transition(structure::ExplicitGSPN, transition_name, transition,
        stoichiometry, dependencies)
    if structure.transition_base<0
        structure.transition_base=length(structure.id_to_pt)
    end
    push!(structure.id_to_pt, transition_name)
    id=length(structure.id_to_pt)
    add_vertex!(structure.gspn, id)
    push!(structure.transition_prop, transition)

    entry_idx=1
    for entry in stoichiometry
        (place, stoichiometric_number)=entry[1:2]
        properties={"stoich"=>stoichiometric_number}
        local_name=string(entry_idx)
        if Base.length(entry)>2
            local_name=entry[3]
        end
        place_id=structure.pt_to_id[place]
        ep=EdgeProperty(stoichiometric_number, local_name)
        push!(structure.edge_prop, ep)
        add_edge!(structure.gspn, id, place_id)
        entry_idx+=1
    end
    # The dependency graph is a separate entity but
    # saved in the same structure with stoichiometry=0.
    for (place, local_name) in dependencies
        place_id=structure.pt_to_id[place]
        add_edge!(structure.gspn, id, place_id)
        push!(structure.edge_prop, EdgeProperty(0, local_name))
    end
end

function transition_places(f::Function, eg::ExplicitGSPN, transition_id::Int64)
    for edge in out_edges(transition_id, eg.gspn)
        ep=eg.edge_prop[edge_index(edge)]
        f(target(edge, eg.gspn), ep.stoichiometry, ep.local_name)
    end
end

function dependent_transitions(eg::ExplicitGSPN, place_id::Set{Int64})
    affected=Set{Int64}()
    for p in place_id
        for on in out_neighbors(p, eg.gspn)
            push!(affected, on)
        end
    end
    affected
end

abstract Model

# A model has both state and structure.
# This is what is presented to the propagator.
type ExplicitGSPNModel <: Model
    structure::ExplicitGSPN
    state
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
    # Policy: Whether the output places must have no tokens.
    transition_places(model.structure, transition_id) do place, stoich, local_name
        if stoich<0
            token_cnt=length(model.state.marking, place)
            @debug("stoich ", transition_id, " (", stoichiometric_number, ", ", token_cnt, ")")
            if (stoich+token_cnt)<0
                @debug("stoich ", transition_id, "false")
                return false
            end
        end
    end
    return true
end

# The correct enabling time to send is the current time if we are asking whether
# this transition is newly-enabled. If we are just getting the distribution of
# an enabled transition, then send its already-known enabling time.
function transition_distribution(model::ExplicitGSPNModel, transition_id, enabling)
    local_marking=Dict{ASCIIString,Any}()
    transition_places(model.structure, transition_id) do place, stoich, name
        if stoich==0
            #@debug("trans ",place," prop ",edge_properties)
            mp=model.state.marking[place]
            local_marking[name]=mp
        end
    end
    transition=transitionobj(model.structure, transition_id)
    dist, invariant=distribution(transition, local_marking, enabling)
end

function examine_transition(model::ExplicitGSPNModel, transition_id,
        enable::Function, disable::Function)
    stoichiometrically_allowed=stoichiometry_satisfied(model, transition_id)
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
end

init(model::ExplicitGSPNModel)=enable_transitions(model)

# Enables transitions consistent with the current marking.
function enable_transitions(model::ExplicitGSPNModel)
  for transition_id in transitionids(model.structure)
    examine_transition(model, transition_id, (x...)->nothing, (x...)->nothing)
  end
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
    tokens=Dict{ASCIIString,Any}() # Map from the edge to the tokens at a place on the edge.
    transition_places(model.structure, transition_id) do place, stoich, name
        if stoich<0
            tokens[name]=empty_container(model.state.marking)
            take!(model.state.marking, place, tokens[name],
                    -stoich)
            push!(affected_places, place)
        end
    end
    
    @debug("SimpleEGSPN.fire in-tokens ", tokens)
    fire_function=transitionobj(model.structure, transition_id).fire
    fire_function(tokens)
    @debug("SimpleEGSPN.fire out-tokens ", tokens)
    all_tokens=empty_container(model.state.marking)
    for (x,y) in tokens
        move!(y, all_tokens, length(y))
    end
    @debug("SimpleEGSPN.fire all-tokens ", all_tokens)
    transition_places(model.structure, transition_id) do place, stoich, name
        if stoich>0
            # Tokens can be created here.
            fill!(model.state.marking, place, all_tokens, stoich)
            push!(affected_places, place)
        end
    end
    # Policy: What to do with too many tokens.
    affected_places
end

function fire(model::ExplicitGSPNModel, id_time::NRTransition,
        enable::Function, disable::Function, rng)
    @debug("ExplicitGSPNModel.fire enter ", model.state.last_fired)
    transition_id=id_time.key
    affected_places=transition_action(model, transition_id)

    mark_disabled!(model.state, transition_id)
    model.state.last_fired=model.structure.id_to_pt[transition_id]
    current_time!(model, id_time.time)

    examine_transitions=dependent_transitions(model.structure, affected_places)
    @debug("ExplicitGSPNModel.fire affected ", affected_places,
            " exam_trans ", examine_transitions)
    for t in examine_transitions
        examine_transition(model, t, enable, disable)
    end
end

function fire(model::ExplicitGSPNModel, id_time::NRTransition, rng)
    fire(model, id_time, (x...)->nothing, (x...)->nothing, rng)
end

