export run_steps

function run_steps(model, sampling, report, rng)
    running=true
    init(sampling, model, rng)
    trans=NRTransition(model.state.last_fired, current_time(model))
    steps=0
    while running
        trans=next(sampling, model, rng)
        if trans.time!=Inf
            fire(sampling, model, trans, rng)
            running=report(model.state)
        else
            running=false
        end
        steps+=1
    end
end
