using DataStructures


cconstruct(x::TypeVar, y)=x
cconstruct(x::DataType, y)=x{y}

# Given a Vector{V}
container_key{V}(::Type{Vector{V}})=Int
container_value{V}(::Type{Vector{V}})=V
container_construct{V}(::Type{Vector{V}})=Array(V, 0)
container_push{V}(c::Vector{V}, v::V)=begin push!(c, v); length(c) end
container_get{V}(c::Vector{V}, k::Int)=c[k]
container_set{V}(c::Vector{V}, k::Int, v::V)=begin c[k]=v; k end
container_iter{V}(c::Vector{V})=1:length(c)

# Given a Deque
container_key{V}(::Type{Deque{V}})=V
container_value{V}(::Type{Deque{V}})=V
container_construct{V}(::Type{Deque{V}})=Deque{V}()
container_push{V}(c::Deque{V}, v::V)=begin push!(c, v); v end
container_get{V}(c::Deque{V}, k::V)=k
container_set{V}(c::Deque{V}, k::V, v::V)=begin c[k]=v; v end
container_iter{V}(c::Deque{V})=c

# Given a Set
container_key{V}(::Type{Set{V}})=V
container_value{V}(::Type{Set{V}})=V
container_construct{V}(::Type{Set{V}})=Set{V}()
container_push{V}(c::Set{V}, v::V)=begin push!(c, v); v end
container_get{V}(c::Set{V}, k::V)=k
container_set{V}(c::Set{V}, k::V, v::V)=begin c[k]=v; v end
container_iter{V}(c::Set{V})=c

# Given a curried Dict, for instance:
# typealias IntDict{V} Dict{Int,V}
container_key{K,V}(::Type{Dict{K,V}})=K
container_value{K,V}(::Type{Dict{K,V}})=V
container_construct{K,V}(::Type{Dict{K,V}})=Dict{K,V}()
container_push{K,V}(c::Dict{K,V}, k::K, v::V)=begin c[k]=v; k end
container_get{K,V}(c::Dict{K,V}, k::Int)=c[k]
container_set{K,V}(c::Dict{K,V}, k::Int, v::V)=begin c[k]=v; k end
container_iter{K,V}(c::Dict{K,V})=keys(c)

