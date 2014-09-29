
# All for the marking.
# Tokens live in containers. Sometimes the container is an integer.
# Tokens are moved, passed around, in the framework, only within
# containers. The user code can look in those containers.
# TokenCreator makes tokens and destroys them. This would handle
# reusing ids, for instance.
# 
import Base: length, pop!, getindex
export TokenMarking, add_tokens, length
export int_marking, getindex
export TokenState, EmptyToken, IdToken

type TokenContainer{T}
    tokens::Array{T,1}
end

function TokenContainer{T}(tokens::Array{T,1})
    TokenContainer(Array(T,0))
end

length(tc::TokenContainer)=length(tc.tokens)

function move!{T}(a::TokenContainer{T}, b::TokenContainer{T}, n::Int)
    append!(b.tokens, splice!(a.tokens, 1:n))
end

type IntContainer
    tokens::Int
end

length(ic::IntContainer)=ic.tokens

function move!(a::IntContainer, b::IntContainer, n)
    a.tokens-=n
    b.tokens+=n
end

abstract TokenCreator

type DefaultTokenCreator{T} <: TokenCreator
    token_type::DataType
    DefaultTokenCreator()=new(T)
end

function create{T}(creator::DefaultTokenCreator{T}, n::Int)
    TokenContainer{T}( [T() for i in 1:n] )
end

function destroy{T}(creator::DefaultTokenCreator{T}, c::TokenContainer{T})
end

pop!{T}(tc::TokenContainer{T})=pop!(tc.tokens)

type IntTokenCreator <: TokenCreator
end

function create(creator::IntTokenCreator, n::Int)
    IntContainer(n)
end

function destroy(creator::IntTokenCreator, c::IntContainer)
end


type EmptyToken
end

type IdToken
    id::Int
end

# Guarantee unique Ids for each token.
type IdTokenCreator <: TokenCreator
    cnt::Int
    token_type::DataType
    IdTokenCreator()=new(0, IdToken)
end

function create(creator::IdTokenCreator, n::Int)
    s=creator.cnt
    creator.cnt+=n
    TokenContainer{IdToken}( [IdToken(i) for i in s:(s+n)] )
end

function destroy(creator::IdTokenCreator, c::TokenContainer{IdToken})
end

# The marking is a dict which represents no tokens
# by not having a key internally.
type TokenMarking{C, TC}
    dict::Dict{Any,C}
    token_creator::TC
end

function TokenMarking(token_type::DataType)
    creator=DefaultTokenCreator{token_type}()
    TokenMarking(Dict{Any,TokenContainer{token_type}}(), creator)
end

function TokenMarking(creator::TokenCreator)
    TokenMarking(Dict{Any,TokenContainer{creator.token_type}}(), creator)
end

function int_marking()
    TokenMarking{IntContainer,IntTokenCreator}(
            Dict{Any,IntContainer}(), IntTokenCreator())
end

function empty_container(marking::TokenMarking)
    create(marking.token_creator, 0)
end

function length(marking::TokenMarking, place)
    if haskey(marking.dict, place)
        return length(marking.dict[place])
    end
    0
end

# Returns a TokenContainer.
function getindex(marking::TokenMarking, place)
    if haskey(marking.dict, place)
        return marking.dict[place]
    else
        marking.dict[place]=create(marking.token_creator, 0)
        return marking.dict[place]
    end
end

function take!(marking::TokenMarking, place, dest, n::Int)
    move!(marking[place], dest, n)
    if length(marking.dict[place])==0
        pop!(marking.dict, place)
    end
end

function fill!(marking::TokenMarking, place, source, n::Int)
    if !haskey(marking.dict, place)
        marking.dict[place]=create(marking.token_creator, 0)
    end
    # This is where we create tokens, using the creator.
    if length(source)<n
        d=n-length(source)
        move!(create(marking.token_creator, d), marking.dict[place], d)
        n=length(source)
    end
    move!(source, marking.dict[place], n)
end

function add_tokens(marking::TokenMarking, place, n::Int)
    if !haskey(marking.dict, place)
        marking.dict[place]=create(marking.token_creator, n)
    else
        move!(create(marking.token_creator, n), marking.dict[place], n)
    end
end


# The state of the system.
type TokenState
    marking # Marking, by place.
    enabling_time::Dict{Any,Float64} # Enabling time, by transition id.
    current_time::Float64 # Current system time.
    # The last transition fired is part of the state only during
    # the brief time invariants about consistency are broken, namely
    # that we update the marking for a transition but have not yet
    # recalculated which transitions are enabled.
    last_fired
end

function TokenState(marking)
    TokenState(marking, Dict{Any,Float64}(), 0.0, nothing)
end
