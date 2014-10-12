include("metapop_herd.jl")

cnt=50
seed=32
if length(ARGS)>0
    cnt=int(ARGS[1])
end
if length(ARGS)>1
    seed=int(ARGS[2])
end
de_params['c']=cnt
disease_nonexponential['c']=cnt

block_cnt=2
block_length=3

rng=MersenneTwister(seed)
println("total individuals ", block_cnt*block_length*2*cnt)
herd_model(disease_nonexponential, 5, 1, 2,
   "metapop.h5", rng)
println("metapop serial about to start second")
# Profile.clear()
# Profile.init(10^7, 0.001)

@time herd_model(disease_nonexponential, cnt, block_cnt, block_length,
   "metapop.h5", rng)

# Profile.print(open("prof.txt","w"))
# Start
# elapsed time: 393.991072396 seconds (12082794888 bytes allocated, 85.22% gc time
# Used graphs with adjacency list and integers
# elapsed time: 47.84666394 seconds (4182068104 bytes allocated, 53.38% gc time)

