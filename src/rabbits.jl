using Logging
@Logging.configure(level=DEBUG)
include("samplesemimarkov.jl")
include("transitiondistributions.jl")
include("partialprocess.jl")

import Base: getindex, setindex!


"""
The board is complicated by having two ways to refer to a space,
either relative to an individual or absolute, so that we
can write hazards that depend on neighbors or on absolute locations.
"""
type Board
  state::Array{Int,2}
  individuals::Array{Tuple{Int,Int}, 1}
  disease::Array{Int, 1}
  M::Int # size of the board
  N::Int # number of individuals
  direction::Array{Int, 2}
  opposite::Array{Int, 1}
  Board(M, N)=new(zeros(Int, M, M), Array(Tuple{Int,Int}, N),
  		zeros(Int, N), M, N,
  		[ [0, -1] [-1, 0] [0, 1] [1, 0] ], [3, 4, 1, 2])
end

Infected(b::Board)=countnz(b.disease)

"""
Are all individuals present and accounted for? This is just for debugging.
"""
function Consistent(b::Board)
	consistent=true
	for c=1:b.N
		xy=b.individuals[c]
		if b.state[xy[1], xy[2]]!=c
			@warn("error individual $c xy $xy")
			consistent=false
			# assert(b.state[xy[1], xy[2]]==c)
		end
	end

	for i=1:b.M
		for j=1:b.M
			if b.state[i, j]>0 && b.individuals[b.state[i, j]]!=(i, j)
				@warn("error state at $i $j")
				consistent=false
				# assert(b.individuals[b.state[i, j]]==(i, j))
			end
		end
	end
	consistent
end

"""
`getindex` is complicated because we refer to each cell in two
different ways, either relative to an individual or in absolute
space. The relative space is the (individual id, -direction).
The absolute space is (x, y). The values are either
0 for nobody there, -1 for off the board, or the index
of the individual in that location.
c=0 for location. c=1 for disease state
"""
function getindex(board::Board, a, b, c)
	if c==1
		if b<0
			individual=board[a, b, 0]
			assert(individual>0)
			return board.disease[individual]
		else
			return board.disease[a]
		end
	elseif b<0
		# Look at relative locations.
		center=board.individuals[a]
		location=[center[1], center[2]] + board.direction[:, -b]
		if !(0 < location[1] < board.M+1)
			return -1
		elseif !(0 < location[2]< board.M+1)
			return -1
		end
		return board.state[location[1], location[2]]
	elseif b==0
		return board.individuals[a]
	else
		# Or just return the actual location.
		return board.state[a, b]
	end
end


"""
`setindex` takes two kinds of arguments. Either tell it
the value to put at a grid cell between (1,1) and (M, M),
where the value is 0 or an individual index, or make the value
an individual index and the grid square a relative location
of the form (individual index, relative location).
c=0 for location. c=1 for disease state.
"""
function setindex!(board::Board, v::Int, a, b, c)
	if c==1
		if b<0
			# a is the individual, b is the relative location.
			# v is the new disease state
			susceptible=board[a, b, 0]
			assert(susceptible>0)
			board.disease[susceptible]=v
			@debug("setindex! disease individual $susceptible to $v")
			return v
		else
			board.disease[a]=v
			return v
		end
	elseif b<1
		# Setting a value at a relative location moves that individual.
		assert(v==a)
		location=board.individuals[a]
		newlocation=[location[1], location[2]]+board.direction[:, -b]
		# Cannot move off of the board.
		assert(0 < newlocation[1] < board.M+1)
		assert(0 < newlocation[2] < board.M+1)
		# Cannot move on top of somebody.
		assert(board.state[newlocation[1], newlocation[2]]==0)
		board.state[location[1], location[2]]=0
		board.state[newlocation[1], newlocation[2]]=v
		board.individuals[a]=(newlocation[1], newlocation[2])
		return v
	else
		# Setting a value absolutely changes individuals' neighbors.
		@debug("Setindex! absolute called.")
		if v>0
			oldindividual=board.individuals[v]
			if 0<oldindividual[1]<board.M+1
				board.state[oldindividual[1], oldindividual[2]]=0
			end
			board.state[a, b]=v
			board.individuals[v]=(a, b)
		else
			ind_idx=board.state[a, b]
			if ind_idx>0
				board.individuals[a, b]=(0,0)
				board.state[a, b]=0
			end
		end
		return v
	end
end


function WanderDirection(state, m::Tuple{Int,Int,Int})
	assert(Consistent(state))
	i=m[1]
	surroundings=[state[i, -j, 0] for j = 1:4]
	state[m[1], m[2], m[3]]=m[1]
	newsurroundings=[state[i, -j, 0] for j = 1:4]
	# changed=Array(Tuple{Int, Int}}, 0)
	changed=[(i, 0, 0)]
	for j=1:length(surroundings)
		# Modification to individual's surroundings.
		if (surroundings[j]==0)!=(newsurroundings[j]==0)
			push!(changed, (i, -j, 0))
		end
		# Modification to surroundings of individuals nearby.
		if surroundings[j]>0
			push!(changed, (surroundings[j], -state.opposite[j], 0))
			# Moving can change whom you can infect, too.
			push!(changed, (surroundings[j], -state.opposite[j], 1))
			push!(changed, (surroundings[j], -j, 1))
		end
		if newsurroundings[j]>0
			push!(changed, (newsurroundings[j], -state.opposite[j], 0))
			push!(changed, (newsurroundings[j], -state.opposite[j], 1))
			push!(changed, (newsurroundings[j], -j, 1))
		end
	end

	if !Consistent(state)
		@warn("Inconsistent board Wander exit $m")
		assert(false)
	end
	@debug("Wander changed $changed")
	return changed
end



type WanderIntensity <: Intensity
	distribution::TransitionDistribution
	enabled::Bool
	lastlocation::Tuple{Int, Int}
	WanderIntensity(distribution)=new(distribution, false, (0, 0))
end

function Update!(wi::WanderIntensity, time, state,
		me::Tuple{Int, Int, Int}, m::Tuple{Int, Int, Int})
	assert(Consistent(state))
	modified=:Undefined
	enabled_now=false

	enabled_now=(state[m[1], m[2], m[3]]==0)
	@debug("Update! me $me m $m was $(wi.enabled) now $enabled_now")

	if enabled_now != wi.enabled
		if enabled_now
			prevenable=wi.distribution.parameters[3]
			@debug("Update!.wander time $time previous $prevenable")
			wi.distribution.parameters[3]=time
			modified=:Enabled
		else
			modified=:Disabled
		end
		wi.enabled=enabled_now
	else
		modified=:Unmodified
	end
	assert(Consistent(state))
	modified
end



function Infect(state, w::Tuple{Int, Int, Int})
	whom=state[w[1], w[2], 0]
	# The disease state of the neighbor changed,
	# and the neighbor's own disease state changed.
	# and the relative disease state of us for the neighbor.
	changed=[w, (whom, 0, 1), (whom, -state.opposite[-w[2]], 1)]
	state[w[1], w[2], w[3]]=1
	if state[whom, 0, 1]!=1
		@warn("Infect w $w whom $whom")
		assert(state[whom, 0, 1]==1)
	end
	return changed
end


type InfectIntensity <: Intensity
	distribution::TransitionDistribution
	enabled::Bool
	InfectIntensity()=new(TransitionExponential(10.0, 0.0), false)
end

function Update!(ii::InfectIntensity, time, state, i, w, s)
	modified=:Undefined
	infectious=state[i[1], i[2], i[3]]==1
	haveneighbor=state[w[1], w[2], w[3]]>0
	enabled=( infectious && haveneighbor && state[s[1], s[2], s[3]]==0)
	if enabled != ii.enabled
		if enabled
			ii.distribution.enabling_time=time
			modified=:Enabled
		else
			modified=:Disabled
		end
		ii.enabled=enabled
	else
		modified=:Unmodified
	end
	@debug("InfectIntensity.Update! $i $s i $infectious, n $haveneighbor ",
			"e $enabled")
	modified
end


function MakeBoard(M, N, rng)
	# M=10 # board dimension
	# N=10 # number of individuals
	state=Board(M, N)

	for set_idx =1:N
		searching=true
		while searching
			location=[rand(rng, 1:M), rand(rng, 1:M)]
			@debug("individual $set_idx location $location")
			if state[location[1], location[2], 0]==0
				state[location[1], location[2], 0]=set_idx
				searching=false
			end
		end
	end

	for inf_idx=1:5
		# Infect individual
		state[inf_idx, 0, 1]=1
	end

	process=PartialProcess(state)

	for midx = 1:N
		for direction = 1:4
			hazard=WanderIntensity(TransitionWeibull(1, 2, 0))
			AddTransition!(process,
				hazard, ((midx, 0, 0), (midx, -direction, 0),),
				WanderDirection, ((midx, -direction, 0),),
				"m$midx$direction")

			infect=InfectIntensity()
			AddTransition!(process,
				infect, ((midx, 0, 1), (midx, -direction, 0),
						(midx, -direction, 1)),
				Infect, ((midx, -direction, 1),),
				"i$midx$direction")
		end
	end

	Consistent(state)
	(process, state)
end



function Run()
	rng=MersenneTwister(333333)
	process, state=MakeBoard(10, 10, rng)
	Init(process)
	# #sampler=FirstReaction()
	sampler=NextReactionHazards()
	time, clock=Dynamics(process, sampler, rng)
	while !isinf(time) && Infected(state)<state.N
		if startswith(clock.name, "i")
			print("$time $clock\n$(state.state)\n")
		end
		
		time, clock=Dynamics(process, sampler, rng)
	end
	print("time $time\n")
	print(state.state, "\n")
	print(state.disease, "\n")
end

Run()
