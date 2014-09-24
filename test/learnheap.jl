using DataStructures

type Transition
  key
  time::Float64
end

function <(a::Transition, B::Transition)
	a.time<b.time
end

function >(a::Transition, B::Transition)
    a.time>b.time
end

function ==(a::Transition, B::Transition)
    a.time==b.time
end

a=Transition("infect", 0.3)
b=Transition("recover", 0.7)
println(a<b)
println(a>b)
println(a==b)

heap=mutable_binary_minheap(Transition)
println(typeof(heap))
handles=[push!(heap, a), push!(heap, b)]
println(typeof(handles[1]))
println(top(heap))
update!(heap, handles[2], Transition("recover", 0.2))
println(top(heap))
# Can we pop any one we choose?
# pop!(heap, handles[2]) # Nope
