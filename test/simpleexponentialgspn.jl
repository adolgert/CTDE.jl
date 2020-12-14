module SimpleExponentialGSPN

using Distributions
using SmallGraphs

import GSPNModel: Model, all_transitions, fire, modified_transitions
import GSPNModel: current_time
export SEGSPN, GraphExponentialTransition
export all_transitions, modified_transitions, fire, current_time


###########################################
# First cut. this assumes a compartmental model
# with exponential transitions. Stoichiometry
# is in the firing function.
###########################################
struct SEGSPN <: Model
	marking
	hazards
	firing
	current_time
end

function sir()
		m=[99,5,0]
		h=[m->1/(0.4*m[1]*m[2]), m->1/(0.5*m[2])]
		f=[m->begin
	        m[1]-=1
	        m[2]+=1
			end,
			m->begin
			  m[2]-=1
			  m[3]+=1
			end]
	    SEGSPN(m, h, f, 0.0)
end

current_time(model::SEGSPN)=model.current_time

function all_transitions(fdist::Function, model::SEGSPN, rng)
  for transition_id in 1:2
  	scale=model.hazards[transition_id](model.marking)
  	println("id ", transition_id, " scale ", scale)
  	if scale>0
	  	fdist(transition_id, Distributions.Exponential(scale),0)
	end
  end
end


function fire(model::SEGSPN, id_time)
	println("SimpleExponentialGSPN: before", model.marking)
	model.firing[id_time[1]](model.marking)
	println("SimpleExponentialGSPN: after", model.marking)
end


################################################################
# Second cut. Using a graph for the gspn and distribution
# dependencies. Graph includes stoichiometry.
# In this version, the marking is in the Places in the graph.
# There is no modification of tokens. The marking is a count.
################################################################
struct GraphExponentialTransition
	distribution
	fire
end

function GraphExponentialTransition(dist)
	GraphExponentialTransition(dist, 2)
end


struct GraphExponentialGSPN <: Model
	# Graph with stoichiometry
	gspn
	# How the stochastic variable depends upon marking at places.
	dependencies
	last_fired
	current_time
end

current_time(model::GraphExponentialGSPN)=model.current_time

function sir_graph()
	u=UndirectedMultiGraph()
	add_node(u, "s", {"mark"=>99})
	add_node(u, "i", {"mark"=>5})
	add_node(u, "r", {"mark"=>0})
	infect_transition=GraphExponentialTransition(
		lm->Exponential(1/(0.4*lm["s"]*lm["i"])))
	recover_transition=GraphExponentialTransition(
		lm->Exponential(1/(0.5*lm["i"])))
	add_node(u, "infect", {"transition"=>infect_transition, "enable"=>-1})
	add_node(u, "recover", {"transition"=>recover_transition, "enable"=>-1})

	# The local annotation is a thought that tokens need
	# to be modified somehow. Also, different token types
	# would have different stochiometries.
	add_edge(u, "s", "infect", {"stoch"=>-1, "local"=>"s"})
	add_edge(u, "i", "infect", {"stoch"=>-1, "local"=>"i"})
	add_edge(u, "infect", "i", {"stoch"=>2, "local"=>"o"})
	add_edge(u, "i", "recover", {"stoch"=>-1, "local"=>"i"})
	add_edge(u, "recover", "r", {"stoch"=>1, "local"=>"r"})

	deps=UndirectedGraph()
	# Local is the name with which the transition
	# can access this place's marking.
	add_edge(deps, "infect", "s", {"local"=>"s"})
	add_edge(deps, "infect", "i", {"local"=>"i"})
	add_edge(deps, "recover", "i", {"local"=>"i"})
    add_node(deps, "r")

	GraphExponentialGSPN(u, deps, nothing, 0.0)
end

function is_enabled(gspn, node)
	for (target, edge_properties) in gspn.edge[node]
		stoichiometry=edge_properties["stoch"]
		marking=gspn.node[target]["mark"]
		if stoichiometry+marking<0
			return false
		end
	end
	return true
end

function transition_distribution(gspn, dependencies, transition)
	local_marking=Dict{Any,Any}()
	for (place, edge_properties) in dependencies.edge[transition]
		println("trans ",place," prop ",edge_properties)
		local_marking[edge_properties["local"]]=gspn.node[place]["mark"]
	end
	gspn.node[transition]["transition"].distribution(local_marking)
end

function all_transitions(fdist::Function, model::GraphExponentialGSPN, rng)
  println("SimpleExponentialGSPN.all_transitions enter")
  for (node, dict) in model.gspn.node
  	if haskey(dict, "transition")
  		println("transition ", node, " ", dict)
  		if is_enabled(model.gspn, node)
  			distribution=transition_distribution(model.gspn, model.dependencies,
                node)
  			fdist(node, distribution, 0, rng)
            if model.gspn.node[node]["enable"]<0
                model.gspn.node[node]["enable"]=current_time(model)
            end
        else
            if model.gspn.node[node]["enable"]>=0
                model.gspn.node[node]["enable"]=-1
            end
  		end
  	end
  end
end

function modified_transitions(model::GraphExponentialGSPN, enable::Function,
	  disable::Function, rng)
    println("SimpleExponentialGSPN.modified_transitions enter")
    if model.last_fired!=nothing
        affected_places=Set([x[1] for x in model.gspn.edge[model.last_fired]])
        examine_transitions=Set()
        for p in affected_places
            union!(examine_transitions,Set([t[1] for t in model.gspn.edge[p]]))
        end
        newly_enabled=Set()
        for t in examine_transitions
            was_enabled=model.gspn.node[t]["enable"]>=0
            now_enabled=is_enabled(model.gspn, t)
            if !was_enabled && now_enabled
                model.gspn.node[t]["enable"]=current_time(model)
                dist=transition_distribution(model.gspn, model.dependencies, t)
                enable(t, dist, current_time(model), rng)
                push!(newly_enabled, t)
            elseif was_enabled && !now_enabled
                model.gspn.node[t]["enable"]=-1
                disable(t, current_time(model), rng)
            end
        end
        # Those which remain enabled but depend on changed places
        # also need to have their hazards re-evaluated.
        examine_distributions=Set()
        for p in affected_places
            union!(examine_distributions,Set(keys(model.dependencies.node[p])))
        end
        for ht in examine_distributions
            if is_enabled(model.gspn, ht) && model.gspn.node[ht]["enable"]>=0
                if !(ht in newly_enabled)
                     dist=transition_distribution(model.gspn,
                        model.dependencies, ht)
                    enable(ht, dist, current_time(model), rng)
                    # Do we reset the enabling time? A policy choice.
                end
            end
        end
    else
    	all_transitions(enable, model, rng)
    end
end

function fire(model::GraphExponentialGSPN, id_time)
	transition_id=id_time[1]
	for (place, edge_properties) in model.gspn.edge[transition_id]
		println("SimpleEGSPN.fire transition ", transition_id, " place ", place,
            " edge ", edge_properties)
		model.gspn.node[place]["mark"]+=edge_properties["stoch"]
	end
	model.last_fired=transition_id
	model.current_time=id_time[2]
end

end # SimpleExponentialGSPN module
