using DataStructures

import Base: length, setindex!, getindex, append!
export length, setindex!, getindex, append!, iter
export ArrayContainer, DequeContainer, SetContainer, DictContainer

# This is a set of abstractions over common containers designed to
# let us treat some value as a key that leads to the member in the
# container. For Array, that's an int. For Dict, that's the actual key.
# For all other containers, the key is the member itself. Weird, but
# useful for keeping links into graph structures.

abstract Container{K,V}

key_type{K,V}(c::Container{K,V})=K
value_type{K,V}(c::Container{K,V})=V
key_type{K,V}(::Type{Container{K,V}})=K
value_type{K,V}(c::Type{Container{K,V}})=V

cconstruct(x::TypeVar, y)=x
cconstruct(x::DataType, y)=x{y}

#### ARRAY #################

type ArrayContainer{V} <: Container{Int,V}
    v::Array{V,1}
    ArrayContainer()=new(Array(V, 0))
    ArrayContainer(cnt::Integer)=new(Array(V, cnt))
end

ArrayContainer(t::DataType)=ArrayContainer{t}()
ArrayContainer(t::DataType, cnt::Integer)=ArrayContainer{t}(cnt)

key_type{V}(::Type{ArrayContainer{V}})=Int

# Returns the key for this vertex.
function append!{V}(c::ArrayContainer{V}, x::V)
    push!(c.v, x)
    length(c.v)
end

function getindex{V}(c::ArrayContainer{V}, k::Int)
    c.v[k]
end

function setindex!{V}(c::ArrayContainer{V}, k::Int, x::V)
    c.v[k]=x
end

length{V}(c::ArrayContainer{V})=length(c.v)

iter(c::ArrayContainer)=1:length(c)

a=ArrayContainer{Int}()
a=ArrayContainer(Int, 3)

#### DEQUE ####################

type DequeContainer{V} <: Container{V,V}
    v::Deque{V}
    DequeContainer()=new(Deque{V}())
end

key_type{V}(::Type{DequeContainer{V}})=V
value_type{V}(::Type{DequeContainer{V}})=V

function append!{V}(c::DequeContainer{V}, x::V)
    push!(c.v, x)
    x
end

function getindex{V}(c::DequeContainer{V}, k::V)
    k
end

function setindex!{V}(c::DequeContainer{V}, k::V, x::V)
    push!(c.v, x)
    x
end

length{V}(c::DequeContainer{V})=length(c.v)
iter(c::DequeContainer)=c.v

a=DequeContainer{Int}()
v=append!(a, 7)

#### SET ######################

type SetContainer{V} <: Container{V,V}
    v::Set{V}
    DequeContainer()=new(Set{V}())
end

function append!{V}(c::SetContainer{V}, x::V)
    push!(c.v, x)
    x
end

function getindex{V}(c::SetContainer{V}, k::V)
    k
end

function setindex!{V}(c::SetContainer{V}, k::V, x::V)
    push!(c.v, x)
    x
end

length{V}(c::SetContainer{V})=length(c.v)
iter(c::SetContainer)=c.v

#### DICT #####################
type DictContainer{K,V} <: Container{K,V}
    v::Dict{K,V}
    DictContainer()=new(Dict{K,V}())
end

function append!{K,V}(c::DictContainer{K,V}, k::K, x::V)
    c.v[k]=x
    k
end

function getindex{K,V}(c::DictContainer{K,V}, k::K)
    c.v[k]
end

function setindex{K,V}(c::DictContainer{K,V}, k::K, x::V)
    c.v[k]=x
end

iter(c::DictContainer)=keys(c.v)
