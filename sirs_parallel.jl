include("sirs.jl")


function pmap_task(f, lst)
    np = nprocs()  # determine the number of processes available
    i = start(lst)
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    nextidx() = begin
        if done(lst, i)
            return i
        end
        (i, idx)=next(lst, i)
        idx
    end
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if done(lst, idx)
                            break
                        end
                        produce(remotecall_fetch(p, f, idx))
                    end
                end
            end
        end
    end
    results
end

function herd_single(params, cnt, obs_times, rng)
    contact=complete_contact_graph(cnt)
    model=individual_exponential_graph(params, contact)
    println("obs_times ", obs_times)

    sampling=NextReactionHazards()
    observer=HerdDiseaseObserver(cnt, obs_times)
    model.state=TokenState(int_marking())
    initialize_marking(model, contact)
    run_to_stop(model, sampling, s->observe(observer, s), rng)
    observer.observations_at_times
end


individual_cnt=int(ARGS[1])
run_cnt=int(ARGS[2])
beta=float(ARGS[3])
gamma=float(ARGS[4])
wane=float(ARGS[5])
seed=int(ARGS[6])
disease_exponential={
    'i'=>beta/individual_cnt, # rate of infection of neighbor
    'r'=>gamma, # infectious to removed
    'w'=>wane, # removed to susceptible
}
observation_times=Time[5.0, 10.0, 15.0]

rng=MersenneTwister(seed)
f(idx)=herd_single(disease_exponential, individual_cnt, observation_times, rng)
f(1)
