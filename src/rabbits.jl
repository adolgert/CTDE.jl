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
  M::Int
  N::Int
  Board(M, N)=new(zeros(Int, M, M), Array(Tuple{Int,Int}, N), M, N)
end

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
"""
function getindex(board::Board, a, b)
	# Look at relative locations.
	if b<0
		direction=Dict(
			 0 => [0, 0],
			-1 => [0, -1],
			-2 => [-1, 0],
			-3 => [0, 1],
			-4 => [1, 0]
			)
		center=board.individuals[a]
		location=[center[1], center[2]] + direction[b]
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
"""
function setindex!(board::Board, v::Int, a, b)
	# Setting a value at a relative location moves that individual.
	if b<1
		assert(v==a)
		direction=Dict(
			 0 => [0, 0],
			-1 => [0, -1],
			-2 => [-1, 0],
			-3 => [0, 1],
			-4 => [1, 0]
			)
		location=board.individuals[a]
		newlocation=[location[1], location[2]]+direction[b]
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


function WanderDirection(state, m::Tuple{Int,Int})
	assert(Consistent(state))
	i=m[1]
	surroundings=[state[i,-1], state[i, -2], state[i, -3], state[i, -4]]
	state[m[1], m[2]]=m[1]
	newsurroundings=[state[i,-1], state[i, -2], state[i, -3], state[i, -4]]
	# changed=Array(Tuple{Int, Int}}, 0)
	changed=[(i, 0)]
	for j=1:length(surroundings)
		# Modification to individual's surroundings.
		if (surroundings[j]==0)!=(newsurroundings[j]==0)
			push!(changed, (i, -j))
		end
		# Modification to surroundings of individuals nearby.
		if surroundings[j]>0
			if j==1
				push!(changed, (surroundings[j], -3))
			elseif j==2
				push!(changed, (surroundings[j], -4))
			elseif j==3
				push!(changed, (surroundings[j], -1))
			elseif j==4
				push!(changed, (surroundings[j], -2))
			end
		end
		if newsurroundings[j]>0
			if j==1
				push!(changed, (newsurroundings[j], -3))
			elseif j==2
				push!(changed, (newsurroundings[j], -4))
			elseif j==3
				push!(changed, (newsurroundings[j], -1))
			elseif j==4
				push!(changed, (newsurroundings[j], -2))
			end
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
		me::Tuple{Int, Int}, m::Tuple{Int, Int})
	assert(Consistent(state))
	modified=:Undefined
	enabled_now=false

	enabled_now=(state[m[1], m[2]]==0)
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



function MakeBoard(M, N, rng)
	# M=10 # board dimension
	# N=10 # number of individuals
	state=Board(M, N)

	for set_idx =1:N
		searching=true
		while searching
			location=[rand(rng, 1:M), rand(rng, 1:M)]
			@debug("individual $set_idx location $location")
			if state[location[1], location[2]]==0
				state[location[1], location[2]]=set_idx
				searching=false
			end
		end
	end

	process=PartialProcess(state)

	for midx = 1:N
		for direction = 1:4
			hazard=WanderIntensity(TransitionWeibull(1, 2, 0))
			action=WanderDirection
			AddTransition!(process, hazard, ((midx, 0), (midx, -direction),),
				WanderDirection, ((midx, -direction),), "$midx$direction")
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
	while !isinf(time) && time<10
		print("$time $clock\n$(state.state)\n")
		
		time, clock=Dynamics(process, sampler, rng)
	end
end

Run()
