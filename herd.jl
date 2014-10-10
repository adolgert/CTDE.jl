include("semimarkov.jl")
using DataFrames
using Gadfly
using SmoothingKernels
using SemiMarkov
import SemiMarkov: enabled_transitions, current_time, current_time!
import SemiMarkov: fire, init


function individual_exponential(params, transition_distributions)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, 'l') # latent
    add_place(structure, 'i') # infectious
    add_place(structure, 'r') # recovered
    add_place(structure, 'n') # not clinical
    add_place(structure, 'c') # clinical

    # The basic SIR
    infectious=ConstExplicitTransition(
        (lm, when)->begin
            (transition_distributions["infectious"](lm, when, params), Int[])
        end)
    recover=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionExponential(params['i'], when), Int[])
        end)
    add_transition(structure, 'f', infectious,
        [('l',-1),('i',1)],
        [])
    add_transition(structure, 'd', recover,
        [('i', -1), ('r', 1)],
        [])

    # The clinical state.
    clinical=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(params['s'], when), Int[])
    	end)
    endclinical=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(params['e'], when), Int[])
    	end)
    add_transition(structure, 'a', clinical,
    	[('n', -1), ('i', -1), ('c', 1), ('i', 1)],
    	[])
    add_transition(structure, 'e', endclinical,
    	[('c', -1), ('r',-1) ('n', 1), ('r',1)],
    	[])

    reset=ConstExplicitTransition( (lm, when)->begin
    	(TransitionExponential(1.0, when), Int[])
	end)
    add_transition(structure, 'z', reset,
    	[('n', -1), ('r', -1), ('n', 1), ('l', 1)],
    	[])


    add_tokens(model.state.marking, 'l', 1)
    add_tokens(model.state.marking, 'n', 1)

    model
end


function individual_nonexponential(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, 'l') # latent
    add_place(structure, 'i') # infectious
    add_place(structure, 'r') # recovered
    add_place(structure, 'n') # not clinical
    add_place(structure, 'c') # clinical

    # The basic SIR
    infectious=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(params['l'][1], params['l'][2], when), Int[])
        end)
    recover=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionGamma(params['i'][1], params['i'][2], when), Int[])
        end)
    add_transition(structure, 'f', infectious,
        [('l',-1),('i',1)],
        [])
    add_transition(structure, 'd', recover,
        [('i', -1), ('r', 1)],
        [])

    # The clinical state.
    clinical=ConstExplicitTransition( (lm, when)->begin
        (TransitionGamma(params['s'][1], params['s'][2], when), Int[])
        end)
    endclinical=ConstExplicitTransition( (lm, when)->begin
        (TransitionExponential(params['e'], when), Int[])
        end)
    add_transition(structure, 'a', clinical,
        [('n', -1), ('i', -1), ('c', 1), ('i', 1)],
        [])
    add_transition(structure, 'e', endclinical,
        [('c', -1), ('r', -1), ('n', 1), ('r', 1)],
        [])

    reset=ConstExplicitTransition( (lm, when)->begin
        (TransitionExponential(1.0, when), Int[])
    end)
    add_transition(structure, 'z', reset,
        [('n', -1), ('r', -1), ('n', 1), ('l', 1)],
        [])


    add_tokens(model.state.marking, 'l', 1)
    add_tokens(model.state.marking, 'n', 1)

    model
end



function individual_independent(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    structure=model.structure
    add_place(structure, 'i') # infected
    add_place(structure, 'r') # recovered
    add_place(structure, 'n') # not clinical
    add_place(structure, 'c') # clinical
    add_place(structure, 'm') # not infectious
    add_place(structure, 'o') # infectious

    # The basic SIR
    infectious=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionWeibull(params['l'][1], params['l'][2], when), Int[])
        end)
    recover=ConstExplicitTransition(
        (lm, when)->begin
            (TransitionGamma(params['i'][1], params['i'][2], when), Int[])
        end)
    add_transition(structure, 'f', infectious,
        [('i',-1),('m',-1),('i',1),('o',1)],
        [])
    add_transition(structure, 'd', recover,
        [('i', -1), ('o', -1),('r',1),('m',1)],
        [])

    # Incubation
    incubate=ConstExplicitTransition( (lm, when)->begin
        (TransitionLogLogistic(params['u'][1], params['u'][2], when), Int[])
        end)
    endclinical=ConstExplicitTransition( (lm, when)->begin
        (TransitionExponential(params['e'], when), Int[])
        end)
    add_transition(structure, 'a', incubate,
        [('n', -1), ('i', -1), ('c', 1), ('i', 1)],
        [])
    add_transition(structure, 'e', endclinical,
        [('c', -1), ('r',-1), ('n', 1), ('r',1)],
        [])

    reset=ConstExplicitTransition( (lm, when)->begin
        (TransitionExponential(1.0, when), Int[])
    end)
    add_transition(structure, 'z', reset,
        [('n', -1), ('r', -1), ('n', 1), ('l', 1)],
        [])

    add_tokens(model.state.marking, 'i', 1)
    add_tokens(model.state.marking, 'n', 1)
    add_tokens(model.state.marking, 'm', 1)

    model
end



function individual_exponential_graph(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    cnt=params['c']
    structure=model.structure
    for i in 1:cnt
        add_place(structure, (i,'s')) # susceptible
        add_place(structure, (i,'l')) # latent
        add_place(structure, (i,'i')) # infectious
        add_place(structure, (i,'r')) # recovered
        add_place(structure, (i,'n')) # not clinical
        add_place(structure, (i,'c')) # clinical
    end

    # The basic SIR
    for i in 1:cnt
        infectious=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionExponential(params['l'], when), Int[])
            end)
        recover=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionExponential(params['i'], when), Int[])
            end)
        add_transition(structure, (i,'f'), infectious,
            [((i,'l'),-1),((i,'i'),1)],
            [])
        add_transition(structure, (i,'d'), recover,
            [((i,'i'), -1), ((i,'r'), 1)],
            [])

        # The clinical state.
        clinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['s'], when), Int[])
            end)
        endclinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['e'], when), Int[])
            end)
        add_transition(structure, (i,'a'), clinical,
            [((i,'n'), -1), ((i,'i'), -1), ((i,'c'), 1), ((i,'i'), 1)],
            [])
        add_transition(structure, (i,'e'), endclinical,
            [((i,'c'), -1), ((i,'n'), 1)],
            [])
    end

    for i in 1:cnt
        for j in 1:cnt
            infect=ConstExplicitTransition((lm, when)->begin
                (TransitionExponential(params['g']/cnt, when), Int[])
                end)
            add_transition(structure, (i,j,'g'), infect,
            [((i,'i'),-1), ((j,'s'), -1), ((i,'i'),1), ((j,'l'),1)],
            [])
        end
    end

    for i in 1:(cnt-1)
        add_tokens(model.state.marking, (i,'s'), 1)
        add_tokens(model.state.marking, (i,'n'), 1)
    end
    add_tokens(model.state.marking, (cnt,'i'), 1)
    add_tokens(model.state.marking, (cnt,'n'), 1)

    model
end


function individual_nonexponential_graph(params)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    cnt=params['c']
    structure=model.structure
    for i in 1:cnt
        add_place(structure, (i,'s')) # susceptible
        add_place(structure, (i,'l')) # latent
        add_place(structure, (i,'i')) # infectious
        add_place(structure, (i,'r')) # recovered
        add_place(structure, (i,'n')) # not clinical
        add_place(structure, (i,'c')) # clinical
    end

    # The basic SIR
    for i in 1:cnt
        infectious=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionWeibull(params['l'][1], params['l'][2], when), Int[])
            end)
        recover=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionGamma(params['i'][1], params['i'][2], when), Int[])
            end)
        add_transition(structure, (i,i,'f'), infectious,
            [((i,'l'),-1),((i,'i'),1)],
            [])
        add_transition(structure, (i,i,'d'), recover,
            [((i,'i'), -1), ((i,'r'), 1)],
            [])

        # The clinical state.
        clinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionGamma(params['s'][1], params['s'][2], when), Int[])
            end)
        endclinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['e'], when), Int[])
            end)
        add_transition(structure, (i,i,'a'), clinical,
            [((i,'n'), -1), ((i,'i'), -1), ((i,'c'), 1), ((i,'i'), 1)],
            [])
        add_transition(structure, (i,i,'e'), endclinical,
            [((i,'c'), -1), ((i,'n'), 1)],
            [])
    end

    for i in 1:cnt
        for j in 1:cnt
            infect=ConstExplicitTransition((lm, when)->begin
                (TransitionExponential(params['g']/cnt, when), Int[])
                end)
            add_transition(structure, (i,j,'g'), infect,
            [((i,'i'),-1), ((j,'s'), -1), ((i,'i'),1), ((j,'l'),1)],
            [])
        end
    end

    for i in 1:(cnt-1)
        add_tokens(model.state.marking, (i,'s'), 1)
        add_tokens(model.state.marking, (i,'n'), 1)
    end
    add_tokens(model.state.marking, (cnt,'i'), 1)
    add_tokens(model.state.marking, (cnt,'n'), 1)

    model
end


# metapopulation of herds
function explicit_metapopulation(params, pen_contact_graph)
    state=TokenState(int_marking())
    model=ExplicitGSPNModel(state)
    pen_cnt=params['c']
    total=pen_cnt*length(pen_contact_graph)
    structure=model.structure
    for i in 1:total
        add_place(structure, (i,'s')) # susceptible
        add_place(structure, (i,'l')) # latent
        add_place(structure, (i,'i')) # infectious
        add_place(structure, (i,'r')) # recovered
        add_place(structure, (i,'n')) # not clinical
        add_place(structure, (i,'c')) # clinical
    end

    # The basic SIR
    for i in 1:total
        infectious=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionWeibull(params['l'][1], params['l'][2], when), Int[])
            end)
        recover=ConstExplicitTransition(
            (lm, when)->begin
                (TransitionGamma(params['i'][1], params['i'][2], when), Int[])
            end)
        add_transition(structure, (i,i,'f'), infectious,
            [((i,'l'),-1),((i,'i'),1)],
            [])
        add_transition(structure, (i,i,'d'), recover,
            [((i,'i'), -1), ((i,'r'), 1)],
            [])

        # The clinical state.
        clinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionGamma(params['s'][1], params['s'][2], when), Int[])
            end)
        endclinical=ConstExplicitTransition( (lm, when)->begin
            (TransitionExponential(params['e'], when), Int[])
            end)
        add_transition(structure, (i,i,'a'), clinical,
            [((i,'n'), -1), ((i,'i'), -1), ((i,'c'), 1), ((i,'i'), 1)],
            [])
        add_transition(structure, (i,i,'e'), endclinical,
            [((i,'c'), -1), ((i,'n'), 1)],
            [])
    end

    for i in 1:total
        peni=1+div(i-1, pen_cnt)
        for j in 1:total
            penj=1+div(j-1, pen_cnt)
            lambda=params['w']/pen_cnt
            if peni==penj
                lambda=params['g']/pen_cnt
            elseif haskey(pen_contact_graph.edge[peni], penj)
                lambda=params['f']/pen_cnt
            end
            infect=ConstExplicitTransition((lm, when)->begin
                (TransitionExponential(lambda, when), Int[])
                end)
            add_transition(structure, (i,j,'g'), infect,
            [((i,'i'),-1), ((j,'s'), -1), ((i,'i'),1), ((j,'l'),1)],
            [])
        end
    end

    for i in 2:total
        add_tokens(model.state.marking, (i,'s'), 1)
        add_tokens(model.state.marking, (i,'n'), 1)
    end
    add_tokens(model.state.marking, (1,'i'), 1)
    add_tokens(model.state.marking, (1,'n'), 1)

    model
end
