
abstract Model

# A model has both state and structure.
# This is what is presented to the propagator.
type GSPNModel{T,S} <: Model
    structure::T
    state::S
end

function current_time(model::GSPNModel)
    model.state.current_time
end

function current_time!(model::GSPNModel, when::Float64)
    @debug("ExplicitModel.current_time! was ",
        model.state.current_time, " is ", when)
    model.state.current_time=when
end

# This exists to translate between user-defined place names
# and the internal place names. It is a facade.
function add_tokens(model::GSPNModel, place::Any, n::Int64)
    place_id=place_to_key(model.structure, place)
    add_tokens(model.state.marking, place_id, n)
end

function stoichiometry_satisfied(model::GSPNModel, transition_id::Int64)
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
    gspnd::ExplicitGSPN
    transition_id::Int
    marking
end

function get(lm::LocalMarking, index::Int)
    lm.marking[transition_dependency(lm.gspnd, lm.transition_id, index)]
end

length(lm::LocalMarking)=transition_dependency_length(lm.gspnd, lm.transition_id)


# The correct enabling time to send is the current time if we are asking whether
# this transition is newly-enabled. If we are just getting the distribution of
# an enabled transition, then send its already-known enabling time.
function transition_distribution(model::GSPNModel, transition_id, enabling)
    local_marking=LocalMarking(model.structure, transition_id,
        model.state.marking)
    transition=transitionobj(model.structure, transition_id)
    dist, invariant=distribution(transition, local_marking,
        model.state.user, enabling)
end

function examine_transition(model::GSPNModel, transition_id,
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

init(model::GSPNModel)=enable_transitions(model)

# Enables transitions consistent with the current marking.
function enable_transitions(model::GSPNModel)
    for transition_id in transitionids(model.structure)
        examine_transition(model, transition_id, (x...)->nothing, (x...)->nothing)
    end
    @debug("enable_transitions enabled ", model.state.enabling)
end


function enabled_transitions(report::Function, model::GSPNModel)
    for (id, enabling) in model.state.enabling
        # This isn't asking who becomes enabled but who is already enabled.
        # So use the already-given enabling time.
        distribution, invariant=transition_distribution(model, id,
                enabling.time)
        @assert(distribution!=nothing)
        report(id, distribution, current_time(model))
    end
    @debug("enabled_transitions enabled ", model.state.enabling)
end


function transition_action(model::GSPNModel, transition_id)
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
    fire(transitionobj(model.structure, transition_id), tokens,
        model.state.user)
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

function fire(model::GSPNModel, id_time::NRTransition,
        enable::Function, disable::Function, rng)
    @debug("GSPNModel.fire enter ", model.state.last_fired, " idt ",
            id_time)
    transition_id=id_time.key
    affected_places=transition_action(model, transition_id)

    mark_disabled!(model.state, transition_id)
    model.state.last_fired=key_to_place(model.structure, transition_id)
    current_time!(model, id_time.time)

    examine_transitions=dependent_transitions(model.structure, affected_places)
    # @debug("GSPNModel.fire affected ", affected_places,
    #         " exam_trans ", examine_transitions)
    for t in examine_transitions
        examine_transition(model, t, enable, disable)
    end
end

function fire(model::GSPNModel, id_time::NRTransition, rng)
    fire(model, id_time, (x...)->nothing, (x...)->nothing, rng)
end
