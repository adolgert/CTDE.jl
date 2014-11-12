# Implements a Finite State Machine, in the category theory
# sense.
export run_steps

function run_steps(model, sampling, report, rng)
    running=true
    init(sampling, model, rng)
    trans=NRTransition(-1, current_time(model))
    steps=0
    while running
        @debug("run_steps next")
        trans=next(sampling, model, rng)
        if trans.time!=Inf
            @debug("run_steps fire")
            fire(sampling, model, trans, rng)
            @debug("run_steps report")
            running=report(model.state)
        else
            @debug("run_steps running=false")
            running=false
        end
        steps+=1
    end
end
