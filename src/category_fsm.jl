# Implements a Finite State Machine, in the category theory
# sense.
export RunSimulation

function RunSimulation(partial_process, sampler, observer, rng)
    running=true
    Init(partial_process)
    steps=0
    while running
        @debug("run_steps next")
        time, clock=Next(sampler, partial_process, rng)
        if !isinf(time)
            running=Fire!(partial_process, time, clock, rng,
                    Observer(sampler), observer)
        else
            @debug("run_steps running=false")
            running=false
        end
        steps+=1
    end
end
