using Distributions
using DataStructures

export FirstReaction, NextReactionHazards
export choose, NRTransition

##############################################
# Classic First Reaction method
##############################################
type FirstReaction
end

function choose(fr::FirstReaction, system, rng)
	first_reaction(system, rng)
end

function first_reaction(system, rng)
	least=NRTransition(nothing, Inf)
	all_transitions(system, rng) do id, dist, now, randgen
	  trial_time=rand(dist, now, randgen)
	  # println("SampleSemiMarkov.first_reaction trial_time ",
	  # 	trial_time, " id ", id, " now ", now)
	  @assert(trial_time>=now)
	  if trial_time<least.time
	  	least.key=id
	  	least.time=trial_time
	  end
    end
    least
end


##################################################
# Next reaction by Hazards
# Also called Anderson's method.
##################################################

# A heap of these records which event comes next.
type NRTransition
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

# Transitions in this method are long-lived. This records historical state.
type TransitionRecord
	remaining_exponential_interval::Float64
	last_modification_time::Float64
	heap_handle::Int64
	distribution
end

# This is the main struct.
type NextReactionHazards
	firing_queue::MutableBinaryHeap{NRTransition,DataStructures.LessThan}
	transition_state::Dict{Any,TransitionRecord}
end

function NextReactionHazards()
	heap=mutable_binary_minheap(NRTransition)
	@debug("SampleSemiMarkov.NextReactionHazards type ",typeof(heap))
	state=Dict{Any,TransitionRecord}()
	NextReactionHazards(heap, state)
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

# Finds the next one and removes it from queue.
function choose(propagator::NextReactionHazards, system, rng)
	@debug("SampleSemiMarkov.choose enter ", current_time(system))
    print_next_reaction_hazards(propagator)
	modified_transitions(system,
		(key, dist, now, randgen)->enable(propagator, key, dist, now, randgen),
		(key, now)->disable(propagator, key, now),
		rng
		)
	now=current_time(system)
	@debug("SampleSemiMarkov.choose ", now, " after modified ",
			propagator.firing_queue)
	key_time=next(propagator, now, rng)
	@debug("SampleSemiMarkov.choose current")
	print_next_reaction_hazards(propagator)
	if key_time.time!=Inf
		# Not firing transition. Removing it from internal queue.
		fire(propagator, key_time.key, now, rng)
	end
	key_time
end

# Finds the next one without removing it from the queue.
function next(propagator::NextReactionHazards, now, rng)
	const NotFound=NRTransition(nothing, Inf)
	if !isempty(propagator.firing_queue)
		least=top(propagator.firing_queue)
	else
		least=NotFound
	end
	@debug("SampleSemiMarkov.next queue length ",
			length(propagator.firing_queue), " least ", least)
	NRTransition(least.key, least.time)
end

function unit_hazard_interval(rng::MersenneTwister)
	-log(rand(rng))
end

# Enable or modify a hazard.
function enable(propagator::NextReactionHazards, key, distribution, now, rng)
	clock_started=haskey(propagator.transition_state, key)
	if clock_started
		record=propagator.transition_state[key]
		when_fire=implicit_hazard_integral(distribution,
			record.remaining_exponential_interval, now)
		@assert(when_fire>=now)
		if record.heap_handle>=0
			@debug("SampleSemiMarkov.enable keyu ", key, " interval ",
				record.remaining_exponential_interval, " when ", when_fire,
				" dist ", distribution)
			update!(propagator.firing_queue, record.heap_handle,
				NRTransition(key, when_fire))
		else
			record.heap_handle=push!(propagator.firing_queue,
				NRTransition(key, when_fire))
			@debug("SampleSemiMarkov.enable keyp ", key, " interval ",
				record.remaining_exponential_interval, " when ", when_fire,
				" dist ", distribution)
		end
		record.last_modification_time=now
		record.distribution=distribution
	else
		interval=unit_hazard_interval(rng)
		firing_time=implicit_hazard_integral(distribution, interval, now)
		@assert(firing_time>=now)
        handle=push!(propagator.firing_queue, NRTransition(key, firing_time))
        @debug("SampleSemiMarkov.enable Adding key ", key, " interval ",
        	interval, " when ", firing_time, " dist ", distribution)
		record=TransitionRecord(interval, now, handle, distribution)
		propagator.transition_state[key]=record
	end
end

# Remove a transition from the queue because it was disabled.
function disable(propagator::NextReactionHazards, key, now)
	record=propagator.transition_state[key]
	# We store distributions in order to calculate remaining hazard
	# which will happen AFTER the state has changed.
	update!(propagator.firing_queue, record.heap_handle,
		NRTransition(key, -1))
	todelete=pop!(propagator.firing_queue)
	@assert(todelete.key==key && todelete.time==-1)

	# Removing the time penalty is what makes this disabling.
	time_penalty=hazard_integral(record.distribution,
		record.last_modification_time, now)
	record.remaining_exponential_interval-=time_penalty
	@debug("SampleSemiMarkov.disable key ", key, " heap length ",
			length(propagator.firing_queue), " time penalty ",
			time_penalty)
	record.last_modification_time=now
	record.distribution=nothing
	record.heap_handle=-1
end

# Remove a transition from the queue because it fired.
function fire(propagator::NextReactionHazards, key, now, rng)
	record=propagator.transition_state[key]
	update!(propagator.firing_queue, record.heap_handle,
		NRTransition(key, -1))
	queue_length=length(propagator.firing_queue)
	removed=pop!(propagator.firing_queue)
	@assert removed.key==key
	@assert queue_length-length(propagator.firing_queue)==1
	@debug("SampleSemiMarkov.fire key ", key, " heap length ",
			length(propagator.firing_queue))
	# Using the same trick for the firing records that we use
	# with the marking. When something is reset, erase it from
	# the dictionary of values. That's equivalent.
	pop!(propagator.transition_state, key)
end

