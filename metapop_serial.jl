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
herd_model(disease_nonexponential, 10, 1, block_length,
   "metapop.h5", rng)
herd_model(disease_nonexponential, cnt, block_cnt, block_length,
   "metapop.h5", rng)
