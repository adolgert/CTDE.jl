include("metapop_herd.jl")

cnt=20
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
block_length=2

rng=MersenneTwister(seed)

herd_model(disease_nonexponential, cnt, block_cnt, block_length,
   "metapop.h5", rng)
