using Distributions
using DataStructures

import Base: <, >, ==
export FirstReaction, NextReactionHazards, NaiveSampler
export NRTransition, Next, Observer

"""
A record of a transition and the time.
It's sortable by time.
"""
immutable NRTransition
	key
	time::Float64
end

function <(a::NRTransition, b::NRTransition)
	a.time<b.time
end

function >(a::NRTransition, b::NRTransition)
    a.time>b.time
end

function ==(a::NRTransition, b::NRTransition)
    a.time==b.time
end



"""
Classic First Reaction method
"""
type FirstReaction
end

function Next(fr::FirstReaction, system, rng)
	least=NRTransition(nothing, Inf)
	Hazards(system, rng) do clock, now, enabled, rng2
	  trial_time=Sample(clock, now, rng2)
	  @assert(trial_time>=now)
	  if trial_time<least.time
	  	least=NRTransition(clock, trial_time)
	  end
    end
    (least.time, least.key)
end

Observer(fr::FirstReaction)=(hazard, time, updated, rng)->nothing


"""
Next reaction by Hazards
Also called Anderson's method.
"""
type TransitionRecord
	exponential_interval::Float64
	heap_handle::Int64
end

type NextReactionHazards
	firing_queue::MutableBinaryHeap{NRTransition,DataStructures.LessThan}
	transition_state::Dict{Any,TransitionRecord}
	init::Bool
end

"""
Construct a Next Reaction sampler.
"""
function NextReactionHazards()
    heap=mutable_binary_minheap(NRTransition)
    @debug("SampleSemiMarkov.NextReactionHazards type ",typeof(heap))
    state=Dict{Any,TransitionRecord}()
    NextReactionHazards(heap, state, true)
end


# Finds the next one without removing it from the queue.
function Next(propagator::NextReactionHazards, system, rng)
	if propagator.init
		Hazards(system, rng) do clock, now, updated, rng2
			Enable(propagator, clock, now, updated, rng2)
	    end
	    propagator.init=false
	end

	const NotFound=NRTransition(nothing, Inf)
	if !isempty(propagator.firing_queue)
		least=top(propagator.firing_queue)
	else
		least=NotFound
	end
	@debug("SampleSemiMarkov.next queue length ",
			length(propagator.firing_queue), " least ", least)
	(least.time, least.key)
end


"""
Returns an observer of intensities to decide what to
do when they change.
"""
function Observer(propagator::NextReactionHazards)
	function nrobserve(clock, time, updated, rng)
		if updated==:Disabled || updated==:Fired
			Disable(propagator, clock, time, updated, rng)
		else
			Enable(propagator, clock, time, updated, rng)
		end
	end
end



function unit_hazard_interval(rng::MersenneTwister)
	-log(rand(rng))
end

# Enable or modify a hazard.
function Enable(propagator::NextReactionHazards, distribution,
		now, updated, rng)
	key=distribution
	clock_started=haskey(propagator.transition_state, key)
	if clock_started
		record=propagator.transition_state[key]
		when_fire=Putative(distribution, now, record.exponential_interval)

		@assert(when_fire>=now)
		if record.heap_handle>=0
			@debug("SampleSemiMarkov.enable keyu ", key, " interval ",
				record.exponential_interval, " when ", when_fire,
				" dist ", distribution)
			update!(propagator.firing_queue, record.heap_handle,
				NRTransition(key, when_fire))
		else
			record.heap_handle=push!(propagator.firing_queue,
				NRTransition(key, when_fire))
			@debug("SampleSemiMarkov.enable keyp ", key, " interval ",
				record.exponential_interval, " when ", when_fire,
				" dist ", distribution)
		end
	else
		interval=unit_hazard_interval(rng)
		firing_time=Putative(distribution, now, interval)
		@assert(firing_time>=now)
        handle=push!(propagator.firing_queue, NRTransition(key, firing_time))
        @debug("SampleSemiMarkov.enable Adding key ", key, " interval ",
        	interval, " when ", firing_time, " dist ", distribution)
		record=TransitionRecord(interval, handle)
		propagator.transition_state[key]=record
	end
    @debug("SampleSemiMarkov.enable exit")
end


# Remove a transition from the queue because it was disabled.
function Disable(propagator::NextReactionHazards, key, now,
        updated, rng)
	record=propagator.transition_state[key]
	# We store distributions in order to calculate remaining hazard
	# which will happen AFTER the state has changed.
	update!(propagator.firing_queue, record.heap_handle,
		NRTransition(key, -1.))
	todelete=pop!(propagator.firing_queue)
	@assert(todelete.key==key && todelete.time==-1)
    if updated==:Disabled
    	record.heap_handle=-1 # This is the official sign it was disabled.
    elseif updated==:Fired
        # Deleting the key is slower for small, finite systems,
        # but it makes infinite (meaning long-running) systems possible.
        delete!(propagator.transition_state, key)
    else
        assert(updated==:Disabled || updated==:Fired)
    end
end


function print_next_reaction_hazards(propagator::NextReactionHazards)
    @debug("NextReactionHazards.firing_queue")
    for n in propagator.firing_queue.nodes
        @debug("  ", n.value)
    end
    arr=Array(Any,0)
    for x in keys(propagator.transition_state)
        push!(arr, x)
    end
    sort!(arr)
    @debug("NextReactionHazards.Transitions")
    @debug("key  remain  last  hazard te")
    for trans in arr
        rec=propagator.transition_state[trans]
        if rec.distribution!=nothing
	        p=parameters(rec.distribution)
    	    @debug(trans, " ", rec.remaining_exponential_interval, " ",
        	    rec.last_modification_time, " ", p[1], " ", p[2])
    	else
    	    @debug(trans, " ", rec.remaining_exponential_interval, " ",
        	    rec.last_modification_time, " nodist")
    	end
    end
end



"""
This implements an incorrect algorithm which is the most common
intuitive choice. Each trajectory's distribution looks correct,
but the master equation isn't correct.
"""

type NaiveSampler
    firing_queue::MutableBinaryHeap{NRTransition,DataStructures.LessThan}
    # This maps from transition to entry in the firing queue.
    transition_entry::Dict{Any,Int}
    init::Bool
end

"""
Construct a NaiveSampler.
"""
function NaiveSampler()
    heap=mutable_binary_minheap(NRTransition)
    @warn("SampleSemiMarkov.NaiveSampler This sampler is WRONG")
    state=Dict{Any,Int}()
    NaiveSampler(heap, state, true)
end


# Finds the next one without removing it from the queue.
function Next(propagator::NaiveSampler, system, rng)
    if propagator.init
        Hazards(system, rng) do clock, now, updated, rng2
            NaiveObserve(propagator, clock, now, updated, rng2)
        end
        propagator.init=false
    end

    const NotFound=NRTransition(nothing, Inf)
    if !isempty(propagator.firing_queue)
        least=top(propagator.firing_queue)
    else
        least=NotFound
    end
    @debug("SampleSemiMarkov.next queue length ",
            length(propagator.firing_queue), " least ", least)
    (least.time, least.key)
end


"""
Returns an observer of intensities to decide what to
do when they change.
"""
function Observer(propagator::NaiveSampler)
    function nrobserve(clock, time, updated, rng)
        NaiveObserve(propagator, clock, time, updated, rng)
    end
end


function NaiveObserve(propagator::NaiveSampler, clock,
            time, updated, rng)
    key=clock
    if updated==:Fired || updated==:Disabled
        heap_handle=propagator.transition_entry[key]
        # We store distributions in order to calculate remaining hazard
        # which will happen AFTER the state has changed.
        update!(propagator.firing_queue, heap_handle,
            NRTransition(key, -1.))
        todelete=pop!(propagator.firing_queue)
        delete!(propagator.transition_entry, key)

    elseif updated==:Enabled
        when_fire=Sample(clock, time, rng)
        heap_handle=push!(propagator.firing_queue,
                NRTransition(key, when_fire))
        propagator.transition_entry[key]=heap_handle

    elseif updated==:Modified
        when_fire=Sample(clock, time, rng)
        heap_handle=propagator.transition_entry[key]
        update!(propagator.firing_queue, heap_handle,
                NRTransition(key, when_fire))
    else
        assert(false)
    end
end
