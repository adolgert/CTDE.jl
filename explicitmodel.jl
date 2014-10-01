
using Distributions
using .SmallGraphs

import Base: push!, length, pop!, get
export ExplicitGSPNModel, ConstExplicitTransition
export ExplicitGSPN, add_place, add_transition, current_time
export sir_explicit, fire


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
        if stoichiometric_number+token_cnt<0
            return false
        end
    end
    return true
end

function transition_distribution(model::ExplicitGSPNModel, transition_id)
    local_marking=Dict{Any,Any}()
    dependencies=model.structure.dependencies
    for (place, edge_properties) in dependencies.edge[transition_id]
        #@debug("trans ",place," prop ",edge_properties)
        mp=model.state.marking[place]
        local_marking[edge_properties["local"]]=mp
    end
    transition=model.structure.gspn.node[transition_id]["transition"]
    distribution(transition, local_marking, current_time(model))
end


# This does two things.
# 1. Enable transitions based on the current state.
# 2. Report them to the callback function.
function all_transitions(enable::Function, disable::Function,
        model::ExplicitGSPNModel, rng)
  for (node, dict) in model.structure.gspn.node
    if haskey(dict, "transition")
        distribution=nothing
        if stoichiometry_satisfied(model, node)
            distribution=transition_distribution(model, node)
        end
        @debug("ExplicitGSPNModel.all_transitions transition ", node,
                " ", dict, " ", distribution)
        if distribution!=nothing
            enable(node, distribution, current_time(model), rng)
            if !marked_enabled(model.state, node)
                mark_enabled!(model.state, node, current_time(model))
            end
        else
            if marked_enabled(model.state, node)
                disable(node, current_time(model))
                mark_disabled!(model.state, node)
            end
        end
    end
  end
  @debug("all_transitions enabled ", model.state.enabling_time)
end


function modified_transitions(model::ExplicitGSPNModel, enable::Function,
      disable::Function, rng)
    @debug("ExplicitGSPNModel.modified_transitions enter ",
            model.state.last_fired)
    if true #model.state.last_fired==nothing
        all_transitions(enable, disable, model, rng)
        return
    end

    gspn=model.structure.gspn
    dependencies=model.structure.dependencies
    affected_places=Set([x[1] for x in gspn.edge[model.state.last_fired]])
    examine_transitions=Set()
    for p in affected_places
        # Transitions whose stoichiometry might change
        union!(examine_transitions, Set([t[1] for t in gspn.edge[p]]))
        # Transitions whose hazards might change (or drop to zero)
        union!(examine_transitions, keys(dependencies.edge[p]))
    end

    @debug("ExplicitGSPNModel.modified_transitions affected ",
            affected_places, " exam_trans ", examine_transitions)
    for t in examine_transitions
        was_enabled=marked_enabled(model.state, t)
        can_enable=stoichiometry_satisfied(model, t)
        dist=nothing
        if can_enable
            dist=transition_distribution(model, t)
        end
        now_enabled=(dist!=nothing)
        @debug("ExplicitGSPNModel.modified_transitions t ", t,
                " was ", was_enabled, " now ", now_enabled)
        if !was_enabled && now_enabled
            mark_enabled!(model.state, t, current_time(model))
            enable(t, dist, current_time(model), rng)
        elseif was_enabled && now_enabled
            enable(t, dist, current_time(model), rng)
        elseif was_enabled && !now_enabled
            mark_disabled!(model.state, t)
            disable(t, current_time(model))
        end
    end
end


# This changes the marking but doesn't update enabling_times
# in order to track which transitions are now enabled.
function fire(model::ExplicitGSPNModel, id_time::NRTransition)
    transition_id=id_time.key
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
        end
    end
    # Policy: What to do with too many tokens.

    mark_disabled!(model.state, transition_id)
    model.state.last_fired=transition_id
    current_time!(model, id_time.time)
end
