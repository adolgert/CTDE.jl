
function sir_explicit(cnt, params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    for p in ["s", "i", "r"]
        add_place(structure, p)
    end

    (β, γ)=params
    infect_transition=ConstExplicitTransition(
        (lm, when)->begin
            s, i=[length(lm[x]) for x in ("s", "i")]
            #println("infect hazard ", hazard, " ", length(lm["s"]), ",", length(lm["i"]))
            (TransitionExponential(β*s*i/cnt, when), Int[s, i])
        end)
    recover_transition=ConstExplicitTransition(
        (lm, when)->begin
            hazard=γ*length(lm["i"])
            #println("recover hazard ", hazard, " ", length(lm["i"]))
            (TransitionExponential(hazard, when), Int[length(lm["i"]),])
        end)
    add_transition(structure, "infect", infect_transition,
        [("s",-1),("i",-1),("i",2)],
        [("s","s"),("i","i")])
    add_transition(structure, "recover", recover_transition,
        [("i", -1), ("r", 1)],
        [("i", "i")])

    add_tokens(model.state.marking, "s", cnt-1)
    add_tokens(model.state.marking, "i", 1)

    model
end


function sir_network(contact)
    model=ExplicitGSPNModel()
    structure=model.structure
    for individual in keys(contact.node)
        for p in ["s", "i", "r"]
            add_place(structure, (individual, p))
        end
    end

    infect_transition=ConstExplicitTransition(
        (lm, when)->begin
            hazard=0.5
            (TransitionExponential(hazard, when), Int[])
        end)
    recover_transition=ConstExplicitTransition(
        (lm, when)->begin
            hazard=0.4
            (TransitionExponential(hazard, when), Int[])
        end)
    for who in keys(contact.node)
        for neighbor in keys(contact.edge[who])
            add_transition(structure, (who,neighbor,"infect"),
                    infect_transition,
                    [((neighbor,"s"),-1),((who, "i"),-1),((who, "i"),2)],
                    [])
        end
        add_transition(structure, "recover", recover_transition,
            [((who,"i"), -1), ((who,"r"), 1)],
            [])
    end

    model
end

sirs_birth_death(cnt)=sirs_birth_death(cnt, [400.0, 0.6, 0, 365/14.0, 1/70.0])

struct SIRSToken
end

# SIRS with birth and death with exponential rates.
function sirs_birth_death(cnt, params)
    state=TokenState(int_marking())
    #state=TokenState(TokenMarking(EmptyToken))
    model=ExplicitGSPNModel(state)
    structure=model.structure
    for p in ["s", "i", "r"]
        add_place(structure, p)
    end

    if length(params)!=5
        error("sirs_birth_death missing params argument")
    end
    (β0, β1, θ, γ, μ)=params
    # This balances birth with death because birth is an absolute hazard.
    α=μ*cnt

    infect_transition=ConstExplicitTransition(
        (m, t)->begin
            (s, i, r)=map(c->length(m[c]), ["s", "i", "r"])
            λ=s*i*β0*(1+β1*cos(2π*(t-θ)))/(s+i+r)
            (TransitionExponential(λ, t), Int[s, i, r])
        end)
    recover_transition=ConstExplicitTransition(
        (m, t)->begin
            (TransitionExponential(γ*length(m["i"]), t), Int[length(m["i"]),])
        end)
    birth_transition=ConstExplicitTransition(
        (m, t)->begin
            (TransitionExponential(α, t), Int[])
        end)
    death_transition=ConstExplicitTransition(
        (m, t)->begin
            (TransitionExponential(μ*length(m["d"]), t), Int[length(m("d")),])
        end)

    add_transition(structure, "infect", infect_transition,
        [("s",-1),("i",-1),("i",2)],
        [("s","s"),("i","i"),("r", "r")])
    add_transition(structure, "recover", recover_transition,
        [("i", -1), ("r", 1)],
        [("i", "i")])
    add_transition(structure, "birth", birth_transition,
        [("i", -1), ("i", 2)],
        [])
    add_transition(structure, "deaths", birth_transition,
        [("s", -1)],
        [("s", "d")])
    add_transition(structure, "deathi", birth_transition,
        [("i", -1)],
        [("s", "d")])
    add_transition(structure, "deathr", birth_transition,
        [("r", -1)],
        [("s", "d")])

    susceptible=("s", int((μ+γ)*cnt/β0))
    infected=("i", int(cnt*(β0-μ-γ)/(β0*(μ+γ))))
    recovered=("r", cnt-(susceptible[2]+infected[2]))

    token_idx=1
    marking=model.state.marking
    for (place, initial_cnt) in [susceptible, infected, recovered]
        add_tokens(model.state.marking, place, initial_cnt)
    end

    model
end
