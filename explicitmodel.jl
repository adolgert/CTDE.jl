
using Distributions
using .SmallGraphs

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


type ExplicitGSPN
    # Graph with stoichiometry
    gspn
    # How the stochastic variable depends upon marking at places.
    dependencies
end

function ExplicitGSPN()
    gspn=UndirectedMultiGraph()
    dependencies=UndirectedGraph()
    ExplicitGSPN(gspn, dependencies)
end

#### Construction
function add_place(structure::ExplicitGSPN, place)
    add_node(structure.gspn, place)
    add_node(structure.dependencies, place)
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
function add_transition(structure::ExplicitGSPN, id, transition,
        stoichiometry, dependencies)
    add_node(structure.gspn, id, {"transition"=>transition})
    entry_idx=1
    for entry in stoichiometry
        (place, stoichiometric_number)=entry[1:2]
        properties={"stoich"=>stoichiometric_number}
        if Base.length(entry)>2
            properties["local"]=entry[3]
        else
            properties["local"]=entry_idx
        end
        if stoichiometric_number>0
            add_edge(structure.gspn, id, place, properties)
        elseif stoichiometric_number<0
            add_edge(structure.gspn, place, id, properties)
        else
            assert(stoichiometric_number!=0)
        end
        entry_idx+=1
    end
    add_node(structure.dependencies, id)
    for (place, local_name) in dependencies
        add_edge(structure.dependencies, id, place, {"local"=>local_name})
    end
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

function current_time(model::ExplicitGSPNModel)
    model.state.current_time
end

function current_time!(model::ExplicitGSPNModel, when::Float64)
    @debug("ExplicitModel.current_time! was ",
        model.state.current_time, " is ", when)
    model.state.current_time=when
end

function stoichiometry_satisfied(model::ExplicitGSPNModel, transition_id)
    # Policy: Whether the output places must have no tokens.
    for (target, edge_properties) in model.structure.gspn.edge[transition_id]
        stoichiometric_number=edge_properties["stoich"]
        token_cnt=length(model.state.marking, target)
        @debug("stoich ", transition_id, " (", stoichiometric_number, ", ", token_cnt, ")")
        if (stoichiometric_number+token_cnt)<0
            @debug("stoich ", transition_id, "false")
            return false
        end
    end
    return true
end

# The correct enabling time to send is the current time if we are asking whether
# this transition is newly-enabled. If we are just getting the distribution of
# an enabled transition, then send its already-known enabling time.
function transition_distribution(model::ExplicitGSPNModel, transition_id, enabling)
    local_marking=Dict{Any,Any}()
    dependencies=model.structure.dependencies
    for (place, edge_properties) in dependencies.edge[transition_id]
        #@debug("trans ",place," prop ",edge_properties)
        mp=model.state.marking[place]
        local_marking[edge_properties["local"]]=mp
    end
    transition=model.structure.gspn.node[transition_id]["transition"]
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
            if (length(invariant)==0) || (enabling_record.invariant!=invariant)
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
  for (node, dict) in model.structure.gspn.node
    if haskey(dict, "transition")
        examine_transition(model, node, (x...)->nothing, (x...)->nothing)
    end
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
    affected_places=Set()
    tokens=Dict() # Map from the edge to the tokens at a place on the edge.
    for (place, edge_properties) in model.structure.gspn.edge[transition_id]
        @debug("SimpleEGSPN.fire transition ", transition_id, " place ", place,
            " edge ", edge_properties)
        stoichiometric_number=edge_properties["stoich"]
        if stoichiometric_number<0
            nickname=edge_properties["local"]
            tokens[nickname]=empty_container(model.state.marking)
            take!(model.state.marking, place, tokens[nickname],
                    -stoichiometric_number)
            push!(affected_places, place)
        end
    end
    
    @debug("SimpleEGSPN.fire in-tokens ", tokens)
    fire_function=model.structure.gspn.node[transition_id]["transition"].fire
    fire_function(tokens)
    @debug("SimpleEGSPN.fire out-tokens ", tokens)
    all_tokens=empty_container(model.state.marking)
    for (x,y) in tokens
        move!(y, all_tokens, length(y))
    end
    @debug("SimpleEGSPN.fire all-tokens ", all_tokens)
    for (place, edge_properties) in model.structure.gspn.edge[transition_id]
        stoichiometric_number=edge_properties["stoich"]
        if stoichiometric_number>0
            nickname=edge_properties["local"]
            # Tokens can be created here.
            fill!(model.state.marking, place, all_tokens, stoichiometric_number)
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
    model.state.last_fired=transition_id
    current_time!(model, id_time.time)

    gspn=model.structure.gspn
    dependencies=model.structure.dependencies
    examine_transitions=Set()
    for p in affected_places
        # Transitions whose stoichiometry might change
        union!(examine_transitions, Set([t[1] for t in gspn.edge[p]]))
        # Transitions whose hazards might change (or drop to zero)
        union!(examine_transitions, keys(dependencies.edge[p]))
    end

    @debug("ExplicitGSPNModel.fire affected ", affected_places,
            " exam_trans ", examine_transitions)
    for t in examine_transitions
        examine_transition(model, t, enable, disable)
    end
end

function fire(model::ExplicitGSPNModel, id_time::NRTransition, rng)
    fire(model, id_time, (x...)->nothing, (x...)->nothing, rng)
end

