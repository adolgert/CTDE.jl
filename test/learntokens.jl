function modify(blah)
    blah["a"][1]+=1
end


marking={1=>[7,2], 2=>[4]}
tokens=Dict()
tokens["a"]=Any[]
push!(tokens["a"], pop!(marking[1]))
println(marking)
println("starting tokens ", tokens)
modify(tokens)
println("modified ", marking)
println("resulting tokens ", tokens)
