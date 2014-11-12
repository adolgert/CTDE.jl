include("semimarkov.jl")
using DataFrames
using Gadfly
using SmoothingKernels
using SemiMarkov
using SemiMarkov.SmallGraphs
import SemiMarkov: enabled_transitions, current_time, current_time!
import SemiMarkov: fire, init

typealias Time Float64

function individual_exponential_graph(params, contact::UndirectedGraph)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    cnt=length(contact)
    structure=model.structure
    for i in 1:cnt
        add_place(structure, (i,'s')) # susceptible
        add_place(structure, (i,'i')) # infectious
        add_place(structure, (i,'r')) # recovered
    end

    # The basic SIR
    for (node, properties) in contact.node
        i=node
        recover=ConstExplicitTransition(
            (lm, when::Time)->begin
                (TransitionExponential(params['r'], when), Int[])
            end)
        add_transition(structure, (i, i, 'd'), recover,
            [((i,'i'), -1), ((i,'r'), 1)],
            [])
    end

    for (source, targets) in contact.edge
        for (target, properties) in targets
            infect=ConstExplicitTransition((lm, when::Time)->begin
                (TransitionExponential(params['i'], when), Int[])
                end)
            (i, j)=(source, target)
            add_transition(structure, (i, j,'g'), infect,
                [((i,'i'),-1), ((j,'s'), -1), ((i,'i'),1), ((j,'i'),1)],
                [])
        end
    end
    model
end

function initialize_marking(model, contact)
    first_infection=1 # Could pick this randomly.
    node_idx=1
    for (node, properties) in contact.node
        if node_idx!=first_infection
            add_tokens(model.state.marking, (node,'s'), 1)
        else
            add_tokens(model.state.marking, (node,'i'), 1)
        end
        node_idx+=1
    end
end

typealias TrajectoryEntry (Int64,Int64,Int64,Time)
type TrajectoryStruct
  s::Int64
  i::Int64
  r::Int64
  t::Time
  TrajectoryStruct()=new(0,0,0,0.)
  TrajectoryStruct(s_,i_,r_,t_)=new(s_,i_,r_,t_)
end

function convert(::Type{TrajectoryEntry}, x::TrajectoryStruct)
    (x.s, x.i, x.r, x.t)
end

type HerdDiseaseObserver
    t::Array{TrajectoryEntry,1}
    cnt::Int
    sir::TrajectoryStruct
    previous_time::Float64
    HerdDiseaseObserver(cnt)=new(Array(TrajectoryEntry, 10_000),
            1, TrajectoryStruct(int(cnt-1), 1, 0, 0.), 0.)
end

function observe(eo::HerdDiseaseObserver, state)
    last_fired=state.last_fired
    delta=state.current_time-eo.previous_time
    shiftstore= xo->begin
        xo.sir.t=state.current_time
        xo.t[xo.cnt]=convert(TrajectoryEntry, xo.sir)
        xo.cnt+=1
    end
    if last_fired[3]=='d'
            eo.sir.i-=1
            eo.sir.r+=1
            shiftstore(eo)
    else
        eo.sir.s-=1
        eo.sir.i+=1
        shiftstore(eo)
    end
    if eo.cnt>length(eo.t)
        new_len=2*eo.cnt
        new_t=Array(TrajectoryEntry, new_len)
        new_t[1:length(eo.t)]=eo.t
        eo.t=new_t
    end
    eo.previous_time=state.current_time
    true # run until you can't run any more.
end

function show(eo::HerdDiseaseObserver)
  for i in 1:(eo.cnt-1)
      println(eo.t[i][1], " ", eo.t[i][2], " ", eo.t[i][3], " ", eo.t[i][4])
  end
end

function epidemic_size(eo::HerdDiseaseObserver)
    eo.t[eo.cnt-1][3]
end

function run_to_stop(model, sampling, report, rng)
    running=true
    init(sampling, model, rng)
    trans=NRTransition(model.state.last_fired, current_time(model))
    steps=0
    while running
        trans=next(sampling, model, rng)
        #println(id_time)
        if trans.time!=Inf
            #println("fire ", id_time)
            fire(sampling, model, trans, rng)
            #println("marking ", model.state.marking)
            running=report(model.state)
        else
            running=false
        end
        steps+=1
    end
    #println("steps ", steps, " end time ", current_time(model))
end

function complete_contact_graph(cnt)
    g=UndirectedGraph()
    for i in 1:cnt
        for j in 1:cnt
            if i!=j
                add_edge(g, i, j)
            end
        end
    end
    g
end


function herd_model(params, cnt, run_cnt, rng)
    contact=complete_contact_graph(cnt)
    model=individual_exponential_graph(params, contact)
    for i in 1:run_cnt
        sampling=NextReactionHazards()
        observer=HerdDiseaseObserver(cnt)
        model.state=TokenState(int_marking())
        initialize_marking(model, contact)
        run_to_stop(model, sampling, s->observe(observer, s), rng)
        println(epidemic_size(observer))
    end
end

individual_cnt=int(ARGS[1])
run_cnt=int(ARGS[2])
beta=float(ARGS[3])
gamma=float(ARGS[4])
seed=int(ARGS[5])
disease_exponential={
    'i'=>beta/individual_cnt, # rate of infection of neighbor
    'r'=>gamma, # latent to infectious
}

rng=MersenneTwister(seed)
herd_model(disease_exponential, individual_cnt, run_cnt, rng)
